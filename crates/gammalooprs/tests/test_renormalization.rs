use gammalooprs::{
    dot,
    graph::{Graph, parse::IntoGraph},
    initialisation::test_initialise,
    model::Model,
    processes::{Amplitude, AmplitudeGraph},
    utils::{GS, load_generic_model},
    uv::{
        ApproximationType, RenormalizationPart, RenormalizationPrescriptionSettings,
        UVOrchestrator, UVgenerationSettings,
        settings::{AlphaLoopSettings, MATADSettings, VakintSettings},
    },
};
use idenso::{
    color::{CS, ColorSimplifier},
    dirac::GammaSimplifier,
    shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip},
};
use spenso::shadowing::symbolica_utils::LogPrint;
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    parse, parse_lit,
};
use symbolica_utils::AtomPrintExt;

fn pole_part_uv_settings() -> UVgenerationSettings {
    let undoing_normalization = parse!(
        "(
     𝑖*(𝜋^((4-2*eps)/2))
  * (exp(-EulerGamma))^(eps)\
  * (exp(-logmUVmu-log_mu_sq))^(eps)\
  )^(-1*(n_loops))",
        default_namespace = "vakint"
    );
    UVgenerationSettings {
        softct: false,
        renormalization_prescription: RenormalizationPrescriptionSettings {
            log_divergent: ApproximationType::PolePart,
            massive_power_divergent: ApproximationType::PolePart,
            massless_power_divergent: ApproximationType::PolePart,
            ..Default::default()
        },
        vakint: VakintSettings {
            normalization: undoing_normalization.to_plain_string(),
            ..Default::default()
        },
        ..Default::default()
    }
}
// GammaLoop and RQFT use opposite ghost-propagator signs (+i versus -i).
// Their ghost-gluon momentum and color conventions differ as well, so the net
// graph sign cannot be inferred from the ghost-edge count or applied globally.
pub fn align_to_rqft(atom: &Atom, model: &Model) -> Atom {
    (model
        .apply_parameter_replacement_rules(
            &model.apply_coupling_replacement_rules(&atom.simplify_color().expand()),
        )
        .replace(parse_lit!(gammalooprs::dim))
        .with(parse_lit!(4))
        .collect_factors()
        .simplify_metrics()
        .simplify_gamma()
        .simplify_color()
        .to_dots()
        .replace(CS.tr)
        .with(Atom::num((1, 2)))
        .replace(CS.nc)
        .with(CS.ca)
        .replace(parse!("UFO::aS"))
        .with(parse!("gs").pow(2) / (Atom::var(Symbol::PI) * 4)))
    .expand_num()
    .collect_factors()
    .collect_num()
    .collect_symbol::<i16>(GS.dim_epsilon)
    .collect_factors()
    // .coefficient_list::<i8>(&[Atom::var(GS.dim_epsilon)])
    // .iter()
    // .fold(Atom::Zero, |a, (e, v)| a + e * v)
}
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
#[test]
fn finite_part_quark_lo() {
    test_initialise().unwrap();
    let g: Vec<Graph> = dot!(digraph d1 {
          overall_factor= "+1"
          projector = "spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(1))))/4/3*(Q(0,spenso::mink(4,1))*spenso::gamma(spenso::bis(4,hedge(1)),spenso::bis(4,hedge(2)),spenso::mink(4,1)))"

          in1 [style=invis];
          in1 -> v1:1
          [particle="d" pin="x:@-left"];

          out1 [style=invis];
          v2:2 -> out1
          [particle="d" pin="x:@+right"];

          v1 -> v2 [particle= "d"];
          v2 -> v1 [particle= "g"];
        }
    )
    .unwrap();

    let mut amp = Amplitude::from_graph_list("bub", g).unwrap();

    let model = load_generic_model("sm");

    let a = amp.graphs[0]
        .renormalization_part(&pole_part_uv_settings())
        .unwrap();

    println!("ren part: {:>}", a.log_print(Some(80)));
    // RQFT quark_lo_0_in.h forest order, H = p1.p1*gs^2*cf (cf = 4/3):
    // F0 (direct) / H = +ep^-1.
    // Sum / H = +ep^-1; native GammaLoop / RQFT = +1. The one-loop Vakint sign
    // and the two quark-gluon vertex phase differences cancel.
    insta::assert_snapshot!(
        align_to_rqft(&a,&model)
        .to_bare_ordered_string(),@"cas(2,cof(3))*dot(P(0,mink(4)),P(0,mink(4)))*gs^2*ε^(-1)"
    );
}

#[test]
fn finite_part_ghost_2loop() {
    test_initialise().unwrap();
    let g: Vec<Graph> = dot!(digraph d1 {//0
      overall_factor= "+1"
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"

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
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"

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
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
      num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
        exte0	 [style=invis];
        exte0	-> 3:0	 [id=0 dir=none particle="ghG" pin="x:@-left"];
        exte1	 [style=invis];
        2:1	-> exte1	 [id=1 dir=none particle="ghG" pin="x:@+right"];
        0:2	-> 1:3	 [id=2 dir=none particle="g"];
        0:4	-> 1:5	 [id=3 dir=none particle="ghG"];
        3:6	-> 0:7	 [id=4 dir=none particle="ghG"];
        1:8	-> 2:9	 [id=5 dir=none particle="ghG"];
        2:10-> 3:11	 [id=6 dir=none particle="g"];
    }

    digraph GL03{//9
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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
        num = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"
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

    let settings = pole_part_uv_settings();

    let new_settings = UVgenerationSettings {
        orchestrator: UVOrchestrator::HedgePoset,
        ..settings.clone()
    };

    let model = load_generic_model("sm");

    #[derive(Debug)]
    struct ForestStatsSnapshot {
        forest_size: usize,
    }

    impl std::fmt::Display for ForestStatsSnapshot {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "forest_size={}", self.forest_size)
        }
    }

    fn assert_new_paths_match_legacy(
        amp: &mut AmplitudeGraph,
        a: RenormalizationPart,
        new_settings: &UVgenerationSettings,
    ) -> ForestStatsSnapshot {
        let normalize = |atom: &Atom| {
            atom.replace(parse_lit!(gammalooprs::dim))
                .with(parse_lit!(4))
                .simplify_metrics()
                .to_dots()
                .simplify_color()
                .expand_num()
                .collect_factors()
        };
        let new_part = amp.renormalization_part(new_settings).unwrap();
        assert_eq!(
            normalize(&new_part.expression),
            normalize(&a.expression),
            "New renormalization for graph {} gives:\n{}\n vs old\n{}",
            amp.graph.name,
            new_part.log_print(Some(120)),
            a.log_print(Some(120))
        );

        ForestStatsSnapshot {
            forest_size: new_part.stats.forest_node_count,
        }
    }

    // Each F_i below follows the summand order in the corresponding RQFT
    // `Fill forest(0)` definition.
    let a = amp.graphs[0].renormalization_part(&settings).unwrap();
    // RQFT ghost_nlo_0_in.h, H = p1.p1*i_*gs^4*ca^2:
    // F0 (140 -> 0) / H = +3/16*ep^-2 + 5/32*ep^-1.
    // F1 (140 -> DM -> 0) / H = -3/16*ep^-2.
    // F2 (140 -> 132 -> 0) / H = -3/16*ep^-2.
    // F3 (not emitted by GammaLoop) / H = 0.
    // Sum / H = -3/16*ep^-2 + 5/32*ep^-1; native GammaLoop / RQFT = +1.
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-3𝑖/16+5𝑖/32*ε)*(cas(2,coad(8)))^2*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-2)");
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[0], a, &new_settings);
    insta::assert_snapshot!(stats.to_string(), @"forest_size=6");

    let a = amp.graphs[1].renormalization_part(&settings).unwrap();
    // RQFT ghost_nlo_1_in.h, H = p1.p1*i_*gs^4*ca^2:
    // F0 (140 -> 0) / H = +1/16*ep^-2 + 1/32*ep^-1.
    // F1 (140 -> oW -> 0) / H = -1/16*ep^-2.
    // F2 (140 -> 132 -> 0) / H = -1/16*ep^-2.
    // F3 (140 -> GS -> 0, not emitted by GammaLoop) / H = 0.
    // Sum / H = -1/16*ep^-2 + 1/32*ep^-1; native GammaLoop / RQFT = -1.
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-1𝑖/32*ε+1𝑖/16)*(cas(2,coad(8)))^2*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-2)"
    );
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[1], a, &new_settings);
    insta::assert_snapshot!(stats.to_string(), @"forest_size=6");

    let a = amp.graphs[2].renormalization_part(&settings).unwrap();
    // RQFT ghost_nlo_2_in.h, H = p1.p1*i_*gs^4*ca^2:
    // F0 (140 -> 0) / H = +1/8*ep^-2 - 1/48*ep^-1.
    // F1 (not emitted by GammaLoop) / H = 0.
    // F2 (140 -> FU -> 0) / H = -1/4*ep^-2 + 1/12*ep^-1.
    // F3 (not emitted by GammaLoop) / H = 0.
    // Sum / H = -1/8*ep^-2 + 1/16*ep^-1; native GammaLoop / RQFT = -1.
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-1𝑖/16*ε+1𝑖/8)*(cas(2,coad(8)))^2*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-2)"
    );
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[2], a, &new_settings);
    insta::assert_snapshot!(stats.to_string(), @"forest_size=4");

    let a = amp.graphs[3].renormalization_part(&settings).unwrap();
    // RQFT ghost_nlo_3_in.h, H = p1.p1*i_*gs^4*ca*nf:
    // F0 (140 -> 0) / H = -1/4*ep^-2 - 13/24*ep^-1.
    // F1 (140 -> zw -> 0) / H = +1/2*ep^-2 + 1/3*ep^-1.
    // Sum / H = +1/4*ep^-2 - 5/24*ep^-1; native GammaLoop / RQFT = +1.
    // The quark-vertex phase cancels the single ghost-propagator convention sign.
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-5𝑖/12*ε+1𝑖/2)*cas(2,coad(8))*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*idx(2,cof(3))*ε^(-2)"
    );
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[3], a, &new_settings);
    insta::assert_snapshot!(stats.to_string(), @"forest_size=4");

    let a = amp.graphs[4].renormalization_part(&settings).unwrap();
    // RQFT ghost_nlo_4_in.h, H = p1.p1*i_*gs^4*ca^2:
    // F0 (140 -> 0) / H = +5/8*ep^-2 - 77/48*ep^-1.
    // F1 (not emitted by GammaLoop) / H = 0.
    // F2 (140 -> zw -> 0) / H = -5/4*ep^-2 + 7/3*ep^-1.
    // F3 (not emitted by GammaLoop) / H = 0.
    // Sum / H = -5/8*ep^-2 + 35/48*ep^-1; native GammaLoop / RQFT = -1.
    // RQFT already includes the graph's 1/2 symmetry factor.
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-35𝑖/48*ε+5𝑖/8)*(cas(2,coad(8)))^2*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-2)"
    );
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[4], a, &new_settings);
    insta::assert_snapshot!(stats.to_string(), @"forest_size=4");

    let a = amp.graphs[5].renormalization_part(&settings).unwrap();
    // RQFT ghost_nlo_5_in.h, H = p1.p1*i_*gs^4*ca^2:
    // F0 (140 -> 0) / H = +5/24*ep^-1.
    // F1 (not emitted by GammaLoop) / H = 0.
    // F2 (140 -> zw -> 0) / H = -1/6*ep^-1.
    // F3 (not emitted by GammaLoop) / H = 0.
    // Sum / H = +1/24*ep^-1; native GammaLoop / RQFT = -1.
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(cas(2,coad(8)))^2*-1𝑖/24*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-1)"
    );
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[5], a, &new_settings);
    insta::assert_snapshot!(stats.to_string(), @"forest_size=4");
}

#[test]
fn finit_part_ghlo() {
    test_initialise().unwrap();

    let g: Vec<Graph> = dot!(digraph d1 {

        projector = "spenso::g(spenso::coad(8,hedge(0)),spenso::coad(8,hedge(1)))/8"

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
        .renormalization_part(&pole_part_uv_settings())
        .unwrap();

    println!("ren part: {:>}", a);
    // RQFT ghost_lo_0_in.h forest order, H = p1.p1*gs^2*ca:
    // F0 (direct) / H = +1/2*ep^-1.
    // Sum / H = +1/2*ep^-1; native GammaLoop / RQFT = +1. GammaLoop's +i
    // ghost propagator differs from RQFT's -i, but the one-loop Vakint sign cancels it.
    insta::assert_snapshot!(
        align_to_rqft(&a,&model)
        .to_bare_ordered_string(),@"1/2*cas(2,coad(8))*dot(P(0,mink(4)),P(0,mink(4)))*gs^2*ε^(-1)"
    );
}

mod failing {
    use super::*;

    fn ghost_3loop_settings() -> UVgenerationSettings {
        UVgenerationSettings {
            softct: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::PolePart,
                massive_power_divergent: ApproximationType::PolePart,
                massless_power_divergent: ApproximationType::PolePart,
                ..Default::default()
            },
            vakint: VakintSettings {
                normalization: "(
                𝑖*(𝜋^((4-2*eps)/2))
             * (exp(-EulerGamma))^(eps)
             * (exp(-logmUVmu-log_mu_sq))^(eps)
             )^(-n_loops)"
                    .to_string(),
                additional_normalization: "1".to_string(),
                matad: MATADSettings {
                    expand_masters: true,
                    susbstitute_masters: true,
                    substitute_hpls: false,
                    direct_numerical_substition: false,
                },
                alphaloop: AlphaLoopSettings {
                    susbstitute_masters: false,
                },
                ..Default::default()
            },
            ..Default::default()
        }
    }

    #[test]
    fn finite_part_ghost_3loop() {
        test_initialise().unwrap();

        let model = load_generic_model("sm");
        let g: Vec<Graph> = Graph::from_path(
            concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/resources/graphs/uv_tests/rqft_ghG_3l.dot"
            ),
            &model,
        )
        .unwrap();

        let mut amp = Amplitude::from_graph_list("bub", g).unwrap();

        let settings = ghost_3loop_settings();

        let a = amp.graphs[0].renormalization_part(&settings).unwrap();
        // `Fi` denotes summand i of RQFT's `Fill forest(0)`.
        // d1: RQFT `ghost_nnlo_0`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(39/16*ep^-2 + 185/32*ep^-1)
        // F1/H = rat(-3/2*ep^-2 - 5/6*ep^-1)
        // F2/H = rat(-3/2*ep^-2 - 5/6*ep^-1)
        // F3/H = rat(-45/16*ep^-2 - 39/8*ep^-1)
        // F4/H = rat(-3/2*ep^-2 - 5/6*ep^-1)
        // F8/H = rat(3/2*ep^-2 + 5/6*ep^-1)
        // F10/H = rat(3/2*ep^-2 + 5/6*ep^-1)
        // F12/H = rat(3/2*ep^-2 + 5/6*ep^-1)
        // F5..F7, F9, F11, F13 = 0
        // sum(Fi)/H = rat(-3/8*ep^-2 + 29/32*ep^-1)
        insta::assert_snapshot!(amp.graphs[0].graph.name,@"d1");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-1456*ε^2+-448/3*ε^2*𝜋^2+-896+4160/3*ε)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        let a = amp.graphs[1].renormalization_part(&settings).unwrap();
        // d2: RQFT `ghost_nnlo_1`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(1/32*ep^-2 + 5/192*ep^-1)
        // F3/H = rat(-3/64*ep^-2)
        // F6/H = rat(-3/64*ep^-2)
        // F1..F2, F4..F5, F7..F16 = 0
        // sum(Fi)/H = rat(-1/16*ep^-2 + 5/192*ep^-1)
        insta::assert_snapshot!(amp.graphs[1].graph.name,@"d2");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"((-32+40/3*ε)*1/64*CA^3+-3/16*f(coad(8,hedge(1)),coad(8,hedge(10)),coad(8,hedge(7)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(8)))*f(coad(8,hedge(10)),coad(8,hedge(12)),coad(8,vertex(4,1)))*f(coad(8,hedge(12)),coad(8,hedge(3)),coad(8,hedge(5)))*f(coad(8,hedge(14)),coad(8,hedge(5)),coad(8,hedge(7)))*f(coad(8,hedge(14)),coad(8,hedge(8)),coad(8,vertex(4,1)))+3/16*f(coad(8,hedge(1)),coad(8,hedge(10)),coad(8,hedge(6)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(8)))*f(coad(8,hedge(10)),coad(8,hedge(12)),coad(8,vertex(4,1)))*f(coad(8,hedge(12)),coad(8,hedge(3)),coad(8,hedge(5)))*f(coad(8,hedge(14)),coad(8,hedge(5)),coad(8,hedge(6)))*f(coad(8,hedge(14)),coad(8,hedge(8)),coad(8,vertex(4,1))))*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-2)"
        );

        let a = amp.graphs[2].renormalization_part(&settings).unwrap();
        // d3: RQFT `ghost_nnlo_2`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(9/128*ep^-3 + 15/256*ep^-2 + 27/512*ep^-1)
        //        + cl2*sqrt3*rat(-3/32*ep^-1)
        //        + pi^2*rat(9/512*ep^-1)
        // F2/H = rat(-27/128*ep^-3 + 9/256*ep^-2)
        //        + pi^2*rat(-9/512*ep^-1)
        // F6/H = rat(-27/128*ep^-3 + 9/256*ep^-2 + 9/512*ep^-1)
        //        + cl2*sqrt3*rat(3/32*ep^-1)
        //        + pi^2*rat(-9/256*ep^-1)
        // F8/H = rat(27/64*ep^-3 - 9/32*ep^-2)
        //        + pi^2*rat(9/256*ep^-1)
        // F1, F3..F5, F7, F9 = 0
        // sum(Fi)/H = rat(9/128*ep^-3 - 39/256*ep^-2 + 9/128*ep^-1)
        insta::assert_snapshot!(amp.graphs[2].graph.name,@"d3");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-52*ε+24+24*ε^2)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        let a = amp.graphs[3].renormalization_part(&settings).unwrap();
        // d4: RQFT `ghost_nnlo_3`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(9/128*ep^-3 + 15/256*ep^-2 + 99/512*ep^-1)
        //        + cl2*sqrt3*rat(-3/32*ep^-1)
        //        + pi^2*rat(9/512*ep^-1)
        // F5/H = rat(-27/128*ep^-3 + 9/256*ep^-2)
        //        + pi^2*rat(-9/512*ep^-1)
        // F7/H = rat(-27/128*ep^-3 + 9/256*ep^-2 + 9/512*ep^-1)
        //        + cl2*sqrt3*rat(3/32*ep^-1)
        //        + pi^2*rat(-9/256*ep^-1)
        // F13/H = rat(27/64*ep^-3 - 9/32*ep^-2)
        //         + pi^2*rat(9/256*ep^-1)
        // F1..F4, F6, F8..F12 = 0
        // sum(Fi)/H = rat(9/128*ep^-3 - 39/256*ep^-2 + 27/128*ep^-1)
        insta::assert_snapshot!(amp.graphs[3].graph.name,@"d4");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-48+56*ε)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-2)"
        );

        let a = amp.graphs[4].renormalization_part(&settings).unwrap();
        // d5: RQFT `ghost_nnlo_4`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(9/32*ep^-2 + 33/64*ep^-1)
        //        + cl2*sqrt3*rat(-3/4*ep^-1)
        // F7/H = rat(-27/64*ep^-2 - 45/128*ep^-1)
        //        + cl2*sqrt3*rat(3/4*ep^-1)
        // F1..F6, F8..F13 = 0
        // sum(Fi)/H = rat(-9/64*ep^-2 + 21/128*ep^-1)
        // The expression snapshot is not yet aligned to this RQFT reference.
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-24+-72*ε^2+52*ε)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        let a = amp.graphs[5].renormalization_part(&settings).unwrap();
        // d6: RQFT `ghost_nnlo_5`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(3/64*ep^-3 - 451/128*ep^-2 - 2767/256*ep^-1)
        //        + cl2*sqrt3*rat(31/16*ep^-1)
        //        + pi^2*rat(3/256*ep^-1)
        // F2/H = rat(27/16*ep^-2 + 19/16*ep^-1)
        // F4/H = rat(-9/64*ep^-3 + 411/128*ep^-2 + 515/64*ep^-1)
        //        + pi^2*rat(-3/256*ep^-1)
        // F5/H = rat(27/16*ep^-2 + 19/16*ep^-1)
        // F8/H = rat(-9/64*ep^-3 + 435/128*ep^-2 + 867/256*ep^-1)
        //        + cl2*sqrt3*rat(-31/16*ep^-1)
        //        + pi^2*rat(-3/128*ep^-1)
        // F10/H = rat(-27/16*ep^-2 - 19/16*ep^-1)
        // F13/H = rat(-27/16*ep^-2 - 19/16*ep^-1)
        // F14/H = rat(9/32*ep^-3 - 45/16*ep^-2 - 13/8*ep^-1)
        //         + pi^2*rat(3/128*ep^-1)
        // F1, F3, F6..F7, F9, F11..F12, F15 = 0
        // sum(Fi)/H = rat(3/64*ep^-3 + 35/128*ep^-2 - ep^-1)
        insta::assert_snapshot!(amp.graphs[5].graph.name,@"d6");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-10712/3*ε+-5576*ε^2+120*ε^2*𝜋^2+3968/3*cl2*sqrt(3)*ε^2+640)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        let a = amp.graphs[6].renormalization_part(&settings).unwrap();
        // d7: RQFT `ghost_nnlo_6`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(3/64*ep^-3 - 451/128*ep^-2 - 2767/256*ep^-1)
        //        + cl2*sqrt3*rat(31/16*ep^-1)
        //        + pi^2*rat(3/256*ep^-1)
        // F2/H = rat(27/16*ep^-2 + 19/16*ep^-1)
        // F4/H = rat(-9/64*ep^-3 + 411/128*ep^-2 + 515/64*ep^-1)
        //        + pi^2*rat(-3/256*ep^-1)
        // F5/H = rat(27/16*ep^-2 + 19/16*ep^-1)
        // F8/H = rat(-9/64*ep^-3 + 435/128*ep^-2 + 867/256*ep^-1)
        //        + cl2*sqrt3*rat(-31/16*ep^-1)
        //        + pi^2*rat(-3/128*ep^-1)
        // F10/H = rat(-27/16*ep^-2 - 19/16*ep^-1)
        // F13/H = rat(-27/16*ep^-2 - 19/16*ep^-1)
        // F14/H = rat(9/32*ep^-3 - 45/16*ep^-2 - 13/8*ep^-1)
        //         + pi^2*rat(3/128*ep^-1)
        // F1, F3, F6..F7, F9, F11..F12, F15 = 0
        // sum(Fi)/H = rat(3/64*ep^-3 + 35/128*ep^-2 - ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-10712/3*ε+-5576*ε^2+-80+3968/3*cl2*sqrt(3)*ε^2)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        let a = amp.graphs[7].renormalization_part(&settings).unwrap();
        // d8: RQFT `ghost_nnlo_7`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(27/32*ep^-3 + 9/32*ep^-2 - 63/64*ep^-1)
        //        + cl2*sqrt3*rat(-9/8*ep^-1)
        //        + pi^2*rat(27/128*ep^-1)
        // F3/H = rat(-81/32*ep^-3 + 27/16*ep^-2)
        //        + pi^2*rat(-27/128*ep^-1)
        // F5/H = rat(-81/64*ep^-3 + 27/128*ep^-2 + 27/256*ep^-1)
        //        + cl2*sqrt3*rat(9/16*ep^-1)
        //        + pi^2*rat(-27/128*ep^-1)
        // F7/H = rat(-81/64*ep^-3 + 27/128*ep^-2 + 27/256*ep^-1)
        //        + cl2*sqrt3*rat(9/16*ep^-1)
        //        + pi^2*rat(-27/128*ep^-1)
        // F10/H = rat(81/32*ep^-3 - 27/16*ep^-2)
        //         + pi^2*rat(27/128*ep^-1)
        // F11/H = rat(81/32*ep^-3 - 27/16*ep^-2)
        //         + pi^2*rat(27/128*ep^-1)
        // F1..F2, F4, F6, F8..F9, F12..F13 = 0
        // sum(Fi)/H = rat(27/32*ep^-3 - 63/64*ep^-2 - 99/128*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-288+264*ε^2+336*ε)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        let a = amp.graphs[8].renormalization_part(&settings).unwrap();
        // d9: RQFT `ghost_nnlo_8`.
        // H = p1.p1*gs^6*ca^3
        // F0..F25 = 0
        // sum(Fi)/H = 0
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(208*z3*ε+60*ε+64)*-1/64*dot(P(0,mink(4)),P(0,mink(4)))*f(coad(8,hedge(1)),coad(8,hedge(15)),coad(8,hedge(9)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(5)))*f(coad(8,hedge(11)),coad(8,hedge(13)),coad(8,hedge(9)))*f(coad(8,hedge(11)),coad(8,hedge(17)),coad(8,hedge(5)))*f(coad(8,hedge(13)),coad(8,hedge(3)),coad(8,hedge(7)))*f(coad(8,hedge(15)),coad(8,hedge(17)),coad(8,hedge(7)))*gs^6*ε^(-2)"
        );

        let a = amp.graphs[9].renormalization_part(&settings).unwrap();
        // d10: RQFT `ghost_nnlo_9`.
        // H = p1.p1*gs^6*ca^3
        // F0..F25 = 0
        // sum(Fi)/H = 0
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"((8*ε+8/3)*1/64*f(coad(8,hedge(1)),coad(8,hedge(11)),coad(8,hedge(15)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(5)))*f(coad(8,hedge(11)),coad(8,hedge(13)),coad(8,hedge(9)))*f(coad(8,hedge(13)),coad(8,hedge(5)),coad(8,hedge(7)))*f(coad(8,hedge(15)),coad(8,hedge(17)),coad(8,hedge(7)))+-1/16*f(coad(8,hedge(1)),coad(8,hedge(11)),coad(8,hedge(15)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(4)))*f(coad(8,hedge(11)),coad(8,hedge(13)),coad(8,hedge(9)))*f(coad(8,hedge(13)),coad(8,hedge(4)),coad(8,hedge(7)))*f(coad(8,hedge(15)),coad(8,hedge(17)),coad(8,hedge(7)))+1/16*f(coad(8,hedge(1)),coad(8,hedge(10)),coad(8,hedge(14)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(5)))*f(coad(8,hedge(10)),coad(8,hedge(13)),coad(8,hedge(9)))*f(coad(8,hedge(13)),coad(8,hedge(5)),coad(8,hedge(7)))*f(coad(8,hedge(14)),coad(8,hedge(17)),coad(8,hedge(7))))*dot(P(0,mink(4)),P(0,mink(4)))*f(coad(8,hedge(17)),coad(8,hedge(3)),coad(8,hedge(9)))*gs^6*ε^(-2)"
        );

        let a = amp.graphs[10].renormalization_part(&settings).unwrap();
        // d11: RQFT `ghost_nnlo_10`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(1/12*ep^-3 + 1/6*ep^-2 + 1/72*ep^-1)
        //        + cl2*sqrt3*rat(-7/27*ep^-1)
        //        + pi^2*rat(1/48*ep^-1)
        // F2/H = rat(-1/8*ep^-3 - 7/48*ep^-2 - 95/864*ep^-1)
        //        + cl2*sqrt3*rat(37/162*ep^-1)
        //        + pi^2*rat(-1/48*ep^-1)
        // F3/H = rat(-1/8*ep^-3 - 5/48*ep^-2 + 35/864*ep^-1)
        //        + cl2*sqrt3*rat(5/162*ep^-1)
        //        + pi^2*rat(-1/48*ep^-1)
        // F6/H = rat(1/4*ep^-3 - 1/36*ep^-1)
        //        + pi^2*rat(1/48*ep^-1)
        // F1, F4..F5, F7 = 0
        // sum(Fi)/H = rat(1/12*ep^-3 - 1/12*ep^-2 - 1/12*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[11].renormalization_part(&settings).unwrap();
        // d12: RQFT `ghost_nnlo_11`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-5/24*ep^-3 + 67/96*ep^-2 + 1049/864*ep^-1)
        //        + cl2*sqrt3*rat(-151/162*ep^-1)
        //        + pi^2*rat(-5/96*ep^-1)
        // F3/H = rat(5/16*ep^-3 - 13/16*ep^-2 - 37/108*ep^-1)
        //        + cl2*sqrt3*rat(185/324*ep^-1)
        //        + pi^2*rat(5/96*ep^-1)
        // F4/H = rat(5/16*ep^-3 - 5/6*ep^-2 - 17/48*ep^-1)
        //        + cl2*sqrt3*rat(13/36*ep^-1)
        //        + pi^2*rat(5/96*ep^-1)
        // F9/H = rat(-5/8*ep^-3 + 115/96*ep^-2 - 95/288*ep^-1)
        //        + pi^2*rat(-5/96*ep^-1)
        // F1..F2, F5..F8, F10..F11 = 0
        // sum(Fi)/H = rat(-5/24*ep^-3 + 1/4*ep^-2 + 3/16*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[12].renormalization_part(&settings).unwrap();
        // d13: RQFT `ghost_nnlo_12`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-7/96*ep^-2 - 113/864*ep^-1)
        //        + cl2*sqrt3*rat(8/81*ep^-1)
        // F3/H = rat(5/96*ep^-2 + 23/1728*ep^-1)
        //        + cl2*sqrt3*rat(-5/162*ep^-1)
        // F4/H = rat(3/32*ep^-2 + 161/1728*ep^-1)
        //        + cl2*sqrt3*rat(-11/162*ep^-1)
        // F9/H = rat(-7/96*ep^-2 + 7/288*ep^-1)
        // F1..F2, F5..F8, F10..F11 = 0
        // sum(Fi)/H = 0
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[13].renormalization_part(&settings).unwrap();
        // d14: RQFT `ghost_nnlo_13`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-3/64*ep^-3 - 29/256*ep^-2 - 335/1536*ep^-1)
        //        + cl2*sqrt3*rat(5/16*ep^-1)
        //        + pi^2*rat(-3/256*ep^-1)
        // F1/H = rat(9/256*ep^-3 - 15/512*ep^-2)
        //        + pi^2*rat(3/1024*ep^-1)
        // F2/H = rat(9/128*ep^-3 + 15/256*ep^-2 + 27/512*ep^-1)
        //        + cl2*sqrt3*rat(-5/32*ep^-1)
        //        + pi^2*rat(3/256*ep^-1)
        // F4/H = rat(9/256*ep^-3 + 21/512*ep^-2)
        //        + pi^2*rat(3/1024*ep^-1)
        // F5/H = rat(9/128*ep^-3 + 15/256*ep^-2 + 27/512*ep^-1)
        //        + cl2*sqrt3*rat(-5/32*ep^-1)
        //        + pi^2*rat(3/256*ep^-1)
        // F9/H = rat(-9/128*ep^-3 + 3/64*ep^-2)
        //        + pi^2*rat(-3/512*ep^-1)
        // F10/H = rat(-9/128*ep^-3 + 3/64*ep^-2)
        //         + pi^2*rat(-3/512*ep^-1)
        // F11/H = rat(-9/128*ep^-3)
        //         + pi^2*rat(-3/512*ep^-1)
        // F3, F6..F8, F12..F14 = 0
        // sum(Fi)/H = rat(-3/64*ep^-3 + 7/64*ep^-2 - 173/1536*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[14].renormalization_part(&settings).unwrap();
        // d15: RQFT `ghost_nnlo_14`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/192*ep^-3 - 1/768*ep^-2 - 5/1536*ep^-1)
        //        + cl2*sqrt3*rat(1/144*ep^-1)
        //        + pi^2*rat(-1/768*ep^-1)
        // F1/H = rat(1/256*ep^-3 - 3/512*ep^-2)
        //        + pi^2*rat(1/3072*ep^-1)
        // F2/H = rat(1/128*ep^-3 + 1/256*ep^-2 + 1/512*ep^-1)
        //        + cl2*sqrt3*rat(-1/288*ep^-1)
        //        + pi^2*rat(1/768*ep^-1)
        // F4/H = rat(1/256*ep^-3 - 3/512*ep^-2)
        //        + pi^2*rat(1/3072*ep^-1)
        // F5/H = rat(1/128*ep^-3 + 1/256*ep^-2 + 1/512*ep^-1)
        //        + cl2*sqrt3*rat(-1/288*ep^-1)
        //        + pi^2*rat(1/768*ep^-1)
        // F9/H = rat(-1/128*ep^-3)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F10/H = rat(-1/128*ep^-3)
        //         + pi^2*rat(-1/1536*ep^-1)
        // F11/H = rat(-1/128*ep^-3)
        //         + pi^2*rat(-1/1536*ep^-1)
        // F3, F6..F8, F12..F14 = 0
        // sum(Fi)/H = rat(-1/192*ep^-3 - 1/192*ep^-2 + 1/1536*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[15].renormalization_part(&settings).unwrap();
        // d16: RQFT `ghost_nnlo_15`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/64*ep^-3 - 13/768*ep^-2 - 53/1536*ep^-1)
        //        + cl2*sqrt3*rat(1/16*ep^-1)
        //        + pi^2*rat(-1/256*ep^-1)
        // F1/H = rat(3/256*ep^-3 - 13/512*ep^-2)
        //        + pi^2*rat(1/1024*ep^-1)
        // F2/H = rat(3/128*ep^-3 + 3/256*ep^-2 + 3/512*ep^-1)
        //        + cl2*sqrt3*rat(-1/96*ep^-1)
        //        + pi^2*rat(1/256*ep^-1)
        // F4/H = rat(3/256*ep^-3 - 1/512*ep^-2)
        //        + pi^2*rat(1/1024*ep^-1)
        // F5/H = rat(3/128*ep^-3 + 5/256*ep^-2 + 9/512*ep^-1)
        //        + cl2*sqrt3*rat(-5/96*ep^-1)
        //        + pi^2*rat(1/256*ep^-1)
        // F9/H = rat(-3/128*ep^-3 + 1/64*ep^-2)
        //        + pi^2*rat(-1/512*ep^-1)
        // F10/H = rat(-3/128*ep^-3)
        //         + pi^2*rat(-1/512*ep^-1)
        // F11/H = rat(-3/128*ep^-3)
        //         + pi^2*rat(-1/512*ep^-1)
        // F3, F6..F8, F12..F14 = 0
        // sum(Fi)/H = rat(-1/64*ep^-3 + 1/384*ep^-2 - 17/1536*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[16].renormalization_part(&settings).unwrap();
        // d17: RQFT `ghost_nnlo_16`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/64*ep^-3 - 19/768*ep^-2 - 21/512*ep^-1)
        //        + cl2*sqrt3*rat(1/16*ep^-1)
        //        + pi^2*rat(-1/256*ep^-1)
        // F1/H = rat(3/256*ep^-3 - 1/512*ep^-2)
        //        + pi^2*rat(1/1024*ep^-1)
        // F2/H = rat(3/128*ep^-3 + 5/256*ep^-2 + 9/512*ep^-1)
        //        + cl2*sqrt3*rat(-5/96*ep^-1)
        //        + pi^2*rat(1/256*ep^-1)
        // F4/H = rat(3/256*ep^-3 - 1/512*ep^-2)
        //        + pi^2*rat(1/1024*ep^-1)
        // F5/H = rat(3/128*ep^-3 + 3/256*ep^-2 + 3/512*ep^-1)
        //        + cl2*sqrt3*rat(-1/96*ep^-1)
        //        + pi^2*rat(1/256*ep^-1)
        // F9/H = rat(-3/128*ep^-3)
        //        + pi^2*rat(-1/512*ep^-1)
        // F10/H = rat(-3/128*ep^-3 + 1/64*ep^-2)
        //         + pi^2*rat(-1/512*ep^-1)
        // F11/H = rat(-3/128*ep^-3)
        //         + pi^2*rat(-1/512*ep^-1)
        // F3, F6..F8, F12..F14 = 0
        // sum(Fi)/H = rat(-1/64*ep^-3 + 7/384*ep^-2 - 9/512*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[17].renormalization_part(&settings).unwrap();
        // d18: RQFT `ghost_nnlo_17`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/128*ep^-3 - 11/768*ep^-2 + 35/512*ep^-1)
        //        + cl2*sqrt3*rat(5/96*ep^-1)
        //        + z3*rat(-1/8*ep^-1)
        //        + pi^2*rat(-1/512*ep^-1)
        // F1/H = rat(3/256*ep^-3 - 1/512*ep^-2)
        //        + pi^2*rat(1/1024*ep^-1)
        // F5/H = rat(3/256*ep^-3 - 13/512*ep^-2)
        //        + pi^2*rat(1/1024*ep^-1)
        // F8/H = rat(3/128*ep^-3 + 5/256*ep^-2 + 9/512*ep^-1)
        //        + cl2*sqrt3*rat(-5/96*ep^-1)
        //        + pi^2*rat(1/256*ep^-1)
        // F11/H = rat(-3/128*ep^-3 + 1/64*ep^-2)
        //         + pi^2*rat(-1/512*ep^-1)
        // F17/H = rat(-3/128*ep^-3 + 1/64*ep^-2)
        //         + pi^2*rat(-1/512*ep^-1)
        // F2..F4, F6..F7, F9..F10, F12..F16, F18 = 0
        // sum(Fi)/H = rat(-1/128*ep^-3 + 7/768*ep^-2 + 11/128*ep^-1)
        //             + z3*rat(-1/8*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[18].renormalization_part(&settings).unwrap();
        // d19: RQFT `ghost_nnlo_18`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-3/128*ep^-3 - 19/256*ep^-2 - 39/512*ep^-1)
        //        + cl2*sqrt3*rat(5/32*ep^-1)
        //        + z3*rat(-1/8*ep^-1)
        //        + pi^2*rat(-3/512*ep^-1)
        // F1/H = rat(9/256*ep^-3 + 21/512*ep^-2)
        //        + pi^2*rat(3/1024*ep^-1)
        // F5/H = rat(9/256*ep^-3 - 15/512*ep^-2)
        //        + pi^2*rat(3/1024*ep^-1)
        // F8/H = rat(9/128*ep^-3 + 15/256*ep^-2 + 27/512*ep^-1)
        //        + cl2*sqrt3*rat(-5/32*ep^-1)
        //        + pi^2*rat(3/256*ep^-1)
        // F11/H = rat(-9/128*ep^-3 + 3/64*ep^-2)
        //         + pi^2*rat(-3/512*ep^-1)
        // F17/H = rat(-9/128*ep^-3 + 3/64*ep^-2)
        //         + pi^2*rat(-3/512*ep^-1)
        // F2..F4, F6..F7, F9..F10, F12..F16, F18 = 0
        // sum(Fi)/H = rat(-3/128*ep^-3 + 23/256*ep^-2 - 3/128*ep^-1)
        //             + z3*rat(-1/8*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[19].renormalization_part(&settings).unwrap();
        // d20: RQFT `ghost_nnlo_19`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/128*ep^-3 - 5/768*ep^-2 - 23/512*ep^-1)
        //        + cl2*sqrt3*rat(1/96*ep^-1)
        //        + pi^2*rat(-1/512*ep^-1)
        // F1/H = rat(3/256*ep^-3 - 1/512*ep^-2)
        //        + pi^2*rat(1/1024*ep^-1)
        // F5/H = rat(3/256*ep^-3 - 1/512*ep^-2)
        //        + pi^2*rat(1/1024*ep^-1)
        // F8/H = rat(3/128*ep^-3 + 3/256*ep^-2 + 3/512*ep^-1)
        //        + cl2*sqrt3*rat(-1/96*ep^-1)
        //        + pi^2*rat(1/256*ep^-1)
        // F11/H = rat(-3/128*ep^-3)
        //         + pi^2*rat(-1/512*ep^-1)
        // F17/H = rat(-3/128*ep^-3)
        //         + pi^2*rat(-1/512*ep^-1)
        // F2..F4, F6..F7, F9..F10, F12..F16, F18 = 0
        // sum(Fi)/H = rat(-1/128*ep^-3 + 1/768*ep^-2 - 5/128*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[20].renormalization_part(&settings).unwrap();
        // d21: RQFT `ghost_nnlo_20`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/384*ep^-3 - 1/256*ep^-2 - 11/1536*ep^-1)
        //        + cl2*sqrt3*rat(1/288*ep^-1)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F1/H = rat(1/256*ep^-3 - 3/512*ep^-2)
        //        + pi^2*rat(1/3072*ep^-1)
        // F5/H = rat(1/256*ep^-3 + 5/512*ep^-2)
        //        + pi^2*rat(1/3072*ep^-1)
        // F8/H = rat(1/128*ep^-3 + 1/256*ep^-2 + 1/512*ep^-1)
        //        + cl2*sqrt3*rat(-1/288*ep^-1)
        //        + pi^2*rat(1/768*ep^-1)
        // F11/H = rat(-1/128*ep^-3)
        //         + pi^2*rat(-1/1536*ep^-1)
        // F17/H = rat(-1/128*ep^-3)
        //         + pi^2*rat(-1/1536*ep^-1)
        // F2..F4, F6..F7, F9..F10, F12..F16, F18 = 0
        // sum(Fi)/H = rat(-1/384*ep^-3 + 1/256*ep^-2 - 1/192*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[21].renormalization_part(&settings).unwrap();
        // d22: RQFT `ghost_nnlo_21`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/128*ep^-3 - 5/768*ep^-2 - 7/512*ep^-1)
        //        + cl2*sqrt3*rat(1/96*ep^-1)
        //        + pi^2*rat(-1/512*ep^-1)
        // F1/H = rat(3/256*ep^-3 - 1/512*ep^-2)
        //        + pi^2*rat(1/1024*ep^-1)
        // F5/H = rat(3/256*ep^-3 - 1/512*ep^-2)
        //        + pi^2*rat(1/1024*ep^-1)
        // F8/H = rat(3/128*ep^-3 + 3/256*ep^-2 + 3/512*ep^-1)
        //        + cl2*sqrt3*rat(-1/96*ep^-1)
        //        + pi^2*rat(1/256*ep^-1)
        // F11/H = rat(-3/128*ep^-3)
        //         + pi^2*rat(-1/512*ep^-1)
        // F17/H = rat(-3/128*ep^-3)
        //         + pi^2*rat(-1/512*ep^-1)
        // F2..F4, F6..F7, F9..F10, F12..F16, F18 = 0
        // sum(Fi)/H = rat(-1/128*ep^-3 + 1/768*ep^-2 - 1/128*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[22].renormalization_part(&settings).unwrap();
        // d23: RQFT `ghost_nnlo_22`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/384*ep^-3 - 1/256*ep^-2 - 11/1536*ep^-1)
        //        + cl2*sqrt3*rat(1/288*ep^-1)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F1/H = rat(1/256*ep^-3 + 5/512*ep^-2)
        //        + pi^2*rat(1/3072*ep^-1)
        // F5/H = rat(1/256*ep^-3 - 3/512*ep^-2)
        //        + pi^2*rat(1/3072*ep^-1)
        // F8/H = rat(1/128*ep^-3 + 1/256*ep^-2 + 1/512*ep^-1)
        //        + cl2*sqrt3*rat(-1/288*ep^-1)
        //        + pi^2*rat(1/768*ep^-1)
        // F11/H = rat(-1/128*ep^-3)
        //         + pi^2*rat(-1/1536*ep^-1)
        // F17/H = rat(-1/128*ep^-3)
        //         + pi^2*rat(-1/1536*ep^-1)
        // F2..F4, F6..F7, F9..F10, F12..F16, F18 = 0
        // sum(Fi)/H = rat(-1/384*ep^-3 + 1/256*ep^-2 - 1/192*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[23].renormalization_part(&settings).unwrap();
        // d24: RQFT `ghost_nnlo_23`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(-1/48*ep^-3 - 5/96*ep^-2 - 7/64*ep^-1)
        //        + cl2*sqrt3*rat(5/36*ep^-1)
        //        + z3*rat(1/4*ep^-1)
        //        + pi^2*rat(-1/192*ep^-1)
        // F1/H = rat(1/32*ep^-3 - 1/64*ep^-2)
        //        + pi^2*rat(1/384*ep^-1)
        // F2/H = rat(1/32*ep^-3 - 1/64*ep^-2)
        //        + pi^2*rat(1/384*ep^-1)
        // F4/H = rat(1/16*ep^-3 + 5/96*ep^-2 + 3/64*ep^-1)
        //        + cl2*sqrt3*rat(-5/36*ep^-1)
        //        + pi^2*rat(1/96*ep^-1)
        // F5/H = rat(-1/16*ep^-3 + 1/24*ep^-2)
        //        + pi^2*rat(-1/192*ep^-1)
        // F6/H = rat(-1/16*ep^-3 + 1/24*ep^-2)
        //        + pi^2*rat(-1/192*ep^-1)
        // F3, F7 = 0
        // sum(Fi)/H = rat(-1/48*ep^-3 + 5/96*ep^-2 - 1/16*ep^-1)
        //             + z3*rat(1/4*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[24].renormalization_part(&settings).unwrap();
        // d25: RQFT `ghost_nnlo_24`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(-1/48*ep^-3 - 5/96*ep^-2 - 7/64*ep^-1)
        //        + cl2*sqrt3*rat(5/36*ep^-1)
        //        + z3*rat(1/4*ep^-1)
        //        + pi^2*rat(-1/192*ep^-1)
        // F1/H = rat(1/32*ep^-3 - 1/64*ep^-2)
        //        + pi^2*rat(1/384*ep^-1)
        // F2/H = rat(1/32*ep^-3 - 1/64*ep^-2)
        //        + pi^2*rat(1/384*ep^-1)
        // F4/H = rat(1/16*ep^-3 + 5/96*ep^-2 + 3/64*ep^-1)
        //        + cl2*sqrt3*rat(-5/36*ep^-1)
        //        + pi^2*rat(1/96*ep^-1)
        // F5/H = rat(-1/16*ep^-3 + 1/24*ep^-2)
        //        + pi^2*rat(-1/192*ep^-1)
        // F6/H = rat(-1/16*ep^-3 + 1/24*ep^-2)
        //        + pi^2*rat(-1/192*ep^-1)
        // F3, F7 = 0
        // sum(Fi)/H = rat(-1/48*ep^-3 + 5/96*ep^-2 - 1/16*ep^-1)
        //             + z3*rat(1/4*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[25].renormalization_part(&settings).unwrap();
        // d26: RQFT `ghost_nnlo_25`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-13/128*ep^-3 - 197/768*ep^-2 - 733/1536*ep^-1)
        //        + cl2*sqrt3*rat(65/96*ep^-1)
        //        + z3*rat(-3/8*ep^-1)
        //        + pi^2*rat(-13/512*ep^-1)
        // F1/H = rat(39/256*ep^-3 - 37/512*ep^-2)
        //        + pi^2*rat(13/1024*ep^-1)
        // F5/H = rat(39/256*ep^-3 - 37/512*ep^-2)
        //        + pi^2*rat(13/1024*ep^-1)
        // F8/H = rat(39/128*ep^-3 + 65/256*ep^-2 + 117/512*ep^-1)
        //        + cl2*sqrt3*rat(-65/96*ep^-1)
        //        + pi^2*rat(13/256*ep^-1)
        // F11/H = rat(-39/128*ep^-3 + 13/64*ep^-2)
        //         + pi^2*rat(-13/512*ep^-1)
        // F17/H = rat(-39/128*ep^-3 + 13/64*ep^-2)
        //         + pi^2*rat(-13/512*ep^-1)
        // F2..F4, F6..F7, F9..F10, F12..F16, F18 = 0
        // sum(Fi)/H = rat(-13/128*ep^-3 + 199/768*ep^-2 - 191/768*ep^-1)
        //             + z3*rat(-3/8*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[26].renormalization_part(&settings).unwrap();
        // d27: RQFT `ghost_nnlo_26`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(1/192*ep^-3 + 7/384*ep^-2 + 31/768*ep^-1)
        //        + cl2*sqrt3*rat(-5/144*ep^-1)
        //        + pi^2*rat(1/768*ep^-1)
        // F1/H = rat(-1/128*ep^-3 - 1/256*ep^-2)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F5/H = rat(-1/128*ep^-3 - 1/256*ep^-2)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F8/H = rat(-1/64*ep^-3 - 5/384*ep^-2 - 3/256*ep^-1)
        //        + cl2*sqrt3*rat(5/144*ep^-1)
        //        + pi^2*rat(-1/384*ep^-1)
        // F11/H = rat(1/64*ep^-3 - 1/96*ep^-2)
        //         + pi^2*rat(1/768*ep^-1)
        // F17/H = rat(1/64*ep^-3 - 1/96*ep^-2)
        //         + pi^2*rat(1/768*ep^-1)
        // F2..F4, F6..F7, F9..F10, F12..F16, F18 = 0
        // sum(Fi)/H = rat(1/192*ep^-3 - 3/128*ep^-2 + 11/384*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[27].renormalization_part(&settings).unwrap();
        // d28: RQFT `ghost_nnlo_27`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/384*ep^-3 - 1/256*ep^-2 - 1/512*ep^-1)
        //        + cl2*sqrt3*rat(5/288*ep^-1)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F1/H = rat(1/256*ep^-3 - 3/512*ep^-2)
        //        + pi^2*rat(1/3072*ep^-1)
        // F5/H = rat(1/256*ep^-3 - 3/512*ep^-2)
        //        + pi^2*rat(1/3072*ep^-1)
        // F8/H = rat(1/128*ep^-3 + 5/768*ep^-2 + 3/512*ep^-1)
        //        + cl2*sqrt3*rat(-5/288*ep^-1)
        //        + pi^2*rat(1/768*ep^-1)
        // F11/H = rat(-1/128*ep^-3 + 1/192*ep^-2)
        //         + pi^2*rat(-1/1536*ep^-1)
        // F17/H = rat(-1/128*ep^-3 + 1/192*ep^-2)
        //         + pi^2*rat(-1/1536*ep^-1)
        // F2..F4, F6..F7, F9..F10, F12..F16, F18 = 0
        // sum(Fi)/H = rat(-1/384*ep^-3 + 1/768*ep^-2 + 1/256*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[28].renormalization_part(&settings).unwrap();
        // d29: RQFT `ghost_nnlo_28`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-3/64*ep^-3 - 7/128*ep^-2 - 23/256*ep^-1)
        //        + cl2*sqrt3*rat(31/144*ep^-1)
        //        + pi^2*rat(-3/256*ep^-1)
        // F1/H = rat(3/64*ep^-3 - 5/128*ep^-2)
        //        + pi^2*rat(1/256*ep^-1)
        // F2/H = rat(3/64*ep^-3 - 1/128*ep^-2 - 11/768*ep^-1)
        //        + cl2*sqrt3*rat(1/144*ep^-1)
        //        + pi^2*rat(1/128*ep^-1)
        // F4/H = rat(3/32*ep^-3 + 3/64*ep^-2 + 25/384*ep^-1)
        //        + cl2*sqrt3*rat(-2/9*ep^-1)
        //        + pi^2*rat(1/64*ep^-1)
        // F7/H = rat(-3/32*ep^-3 + 3/32*ep^-2)
        //        + pi^2*rat(-1/128*ep^-1)
        // F9/H = rat(-3/32*ep^-3 + 1/32*ep^-2)
        //        + pi^2*rat(-1/128*ep^-1)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-3/64*ep^-3 + 9/128*ep^-2 - 5/128*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[29].renormalization_part(&settings).unwrap();
        // d30: RQFT `ghost_nnlo_29`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(1/48*ep^-3 + 1/12*ep^-2 + 25/864*ep^-1)
        //        + cl2*sqrt3*rat(-25/324*ep^-1)
        //        + pi^2*rat(1/192*ep^-1)
        // F1/H = rat(-1/24*ep^-2)
        // F2/H = rat(-1/32*ep^-3 - 13/192*ep^-2 - 23/1152*ep^-1)
        //        + cl2*sqrt3*rat(7/216*ep^-1)
        //        + pi^2*rat(-1/192*ep^-1)
        // F3/H = rat(-1/32*ep^-3 - 5/64*ep^-2 - 103/3456*ep^-1)
        //        + cl2*sqrt3*rat(29/648*ep^-1)
        //        + pi^2*rat(-1/192*ep^-1)
        // F5/H = rat(1/24*ep^-2)
        // F6/H = rat(1/16*ep^-3 + 1/24*ep^-2)
        //        + pi^2*rat(1/192*ep^-1)
        // F4, F7 = 0
        // sum(Fi)/H = rat(1/48*ep^-3 - 1/48*ep^-2 - 1/48*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[30].renormalization_part(&settings).unwrap();
        // d31: RQFT `ghost_nnlo_30`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-19/384*ep^-3 + 131/768*ep^-2 + 4471/13824*ep^-1)
        //        + cl2*sqrt3*rat(-287/2592*ep^-1)
        //        + pi^2*rat(-19/1536*ep^-1)
        // F1/H = rat(-1/128*ep^-3 - 71/768*ep^-2)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F2/H = rat(5/64*ep^-3 - 77/384*ep^-2 - 149/768*ep^-1)
        //        + cl2*sqrt3*rat(11/144*ep^-1)
        //        + pi^2*rat(5/384*ep^-1)
        // F4/H = rat(9/128*ep^-3 - 51/256*ep^-2 - 1105/13824*ep^-1)
        //        + cl2*sqrt3*rat(89/2592*ep^-1)
        //        + pi^2*rat(3/256*ep^-1)
        // F7/H = rat(1/64*ep^-3 + 17/192*ep^-2)
        //        + pi^2*rat(1/768*ep^-1)
        // F9/H = rat(-5/32*ep^-3 + 7/24*ep^-2)
        //        + pi^2*rat(-5/384*ep^-1)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-19/384*ep^-3 + 15/256*ep^-2 + 19/384*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[31].renormalization_part(&settings).unwrap();
        // d32: RQFT `ghost_nnlo_31`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/384*ep^-3 - 19/768*ep^-2 - 719/13824*ep^-1)
        //        + cl2*sqrt3*rat(67/2592*ep^-1)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F1/H = rat(1/128*ep^-3 + 7/768*ep^-2)
        //        + pi^2*rat(1/1536*ep^-1)
        // F2/H = rat(5/192*ep^-2 + 47/1152*ep^-1)
        //        + cl2*sqrt3*rat(-1/54*ep^-1)
        // F4/H = rat(1/128*ep^-3 + 5/256*ep^-2 + 119/13824*ep^-1)
        //        + cl2*sqrt3*rat(-19/2592*ep^-1)
        //        + pi^2*rat(1/768*ep^-1)
        // F7/H = rat(-1/64*ep^-3 - 1/192*ep^-2)
        //        + pi^2*rat(-1/768*ep^-1)
        // F9/H = rat(-1/48*ep^-2)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-1/384*ep^-3 + 1/256*ep^-2 - 1/384*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[32].renormalization_part(&settings).unwrap();
        // d33: RQFT `ghost_nnlo_32`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/64*ep^-3 - 1/128*ep^-2 - 131/6912*ep^-1)
        //        + cl2*sqrt3*rat(25/1296*ep^-1)
        //        + pi^2*rat(-1/256*ep^-1)
        // F1/H = rat(1/64*ep^-3 - 1/384*ep^-2)
        //        + pi^2*rat(1/768*ep^-1)
        // F2/H = rat(1/64*ep^-3 - 1/384*ep^-2 - 11/2304*ep^-1)
        //        + cl2*sqrt3*rat(1/432*ep^-1)
        //        + pi^2*rat(1/384*ep^-1)
        // F4/H = rat(1/32*ep^-3 + 1/192*ep^-2 + 37/3456*ep^-1)
        //        + cl2*sqrt3*rat(-7/324*ep^-1)
        //        + pi^2*rat(1/192*ep^-1)
        // F7/H = rat(-1/32*ep^-3 + 1/96*ep^-2)
        //        + pi^2*rat(-1/384*ep^-1)
        // F9/H = rat(-1/32*ep^-3 + 1/96*ep^-2)
        //        + pi^2*rat(-1/384*ep^-1)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-1/64*ep^-3 + 5/384*ep^-2 - 5/384*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[33].renormalization_part(&settings).unwrap();
        // d34: RQFT `ghost_nnlo_33`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(1/16*ep^-3 + 1/3*ep^-2 + 5/32*ep^-1)
        //        + cl2*sqrt3*rat(-19/36*ep^-1)
        //        + pi^2*rat(1/64*ep^-1)
        // F1/H = rat(-1/8*ep^-2)
        // F2/H = rat(-3/32*ep^-3 - 13/64*ep^-2 - 23/384*ep^-1)
        //        + cl2*sqrt3*rat(7/72*ep^-1)
        //        + pi^2*rat(-1/64*ep^-1)
        // F3/H = rat(-3/32*ep^-3 - 23/64*ep^-2 - 53/384*ep^-1)
        //        + cl2*sqrt3*rat(31/72*ep^-1)
        //        + pi^2*rat(-1/64*ep^-1)
        // F5/H = rat(1/8*ep^-2)
        // F6/H = rat(3/16*ep^-3 + 1/8*ep^-2)
        //        + pi^2*rat(1/64*ep^-1)
        // F4, F7 = 0
        // sum(Fi)/H = rat(1/16*ep^-3 - 5/48*ep^-2 - 1/24*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[34].renormalization_part(&settings).unwrap();
        // d35: RQFT `ghost_nnlo_34`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-19/128*ep^-3 + 229/768*ep^-2 + 721/1536*ep^-1)
        //        + cl2*sqrt3*rat(167/288*ep^-1)
        //        + pi^2*rat(-19/512*ep^-1)
        // F1/H = rat(-3/128*ep^-3 - 67/256*ep^-2)
        //        + pi^2*rat(-1/512*ep^-1)
        // F2/H = rat(15/64*ep^-3 - 77/128*ep^-2 - 149/256*ep^-1)
        //        + cl2*sqrt3*rat(11/48*ep^-1)
        //        + pi^2*rat(5/128*ep^-1)
        // F4/H = rat(27/128*ep^-3 - 69/256*ep^-2 + 217/1536*ep^-1)
        //        + cl2*sqrt3*rat(-233/288*ep^-1)
        //        + pi^2*rat(9/256*ep^-1)
        // F7/H = rat(3/64*ep^-3 + 15/64*ep^-2)
        //        + pi^2*rat(1/256*ep^-1)
        // F9/H = rat(-15/32*ep^-3 + 7/8*ep^-2)
        //        + pi^2*rat(-5/128*ep^-1)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-19/128*ep^-3 + 211/768*ep^-2 + 11/384*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[35].renormalization_part(&settings).unwrap();
        // d36: RQFT `ghost_nnlo_35`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/128*ep^-3 - 53/768*ep^-2 - 83/512*ep^-1)
        //        + cl2*sqrt3*rat(7/96*ep^-1)
        //        + pi^2*rat(-1/512*ep^-1)
        // F1/H = rat(3/128*ep^-3 + 3/256*ep^-2)
        //        + pi^2*rat(1/512*ep^-1)
        // F2/H = rat(5/64*ep^-2 + 47/384*ep^-1)
        //        + cl2*sqrt3*rat(-1/18*ep^-1)
        // F4/H = rat(3/128*ep^-3 + 11/256*ep^-2 + 25/1536*ep^-1)
        //        + cl2*sqrt3*rat(-5/288*ep^-1)
        //        + pi^2*rat(1/256*ep^-1)
        // F7/H = rat(-3/64*ep^-3 + 1/64*ep^-2)
        //        + pi^2*rat(-1/256*ep^-1)
        // F9/H = rat(-1/16*ep^-2)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-1/128*ep^-3 + 13/768*ep^-2 - 3/128*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[36].renormalization_part(&settings).unwrap();
        // d37: RQFT `ghost_nnlo_36`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-3/64*ep^-3 - 7/128*ep^-2 - 23/256*ep^-1)
        //        + cl2*sqrt3*rat(31/144*ep^-1)
        //        + pi^2*rat(-3/256*ep^-1)
        // F1/H = rat(3/64*ep^-3 - 5/128*ep^-2)
        //        + pi^2*rat(1/256*ep^-1)
        // F2/H = rat(3/64*ep^-3 - 1/128*ep^-2 - 11/768*ep^-1)
        //        + cl2*sqrt3*rat(1/144*ep^-1)
        //        + pi^2*rat(1/128*ep^-1)
        // F4/H = rat(3/32*ep^-3 + 3/64*ep^-2 + 25/384*ep^-1)
        //        + cl2*sqrt3*rat(-2/9*ep^-1)
        //        + pi^2*rat(1/64*ep^-1)
        // F7/H = rat(-3/32*ep^-3 + 3/32*ep^-2)
        //        + pi^2*rat(-1/128*ep^-1)
        // F9/H = rat(-3/32*ep^-3 + 1/32*ep^-2)
        //        + pi^2*rat(-1/128*ep^-1)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-3/64*ep^-3 + 9/128*ep^-2 - 5/128*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[37].renormalization_part(&settings).unwrap();
        // d38: RQFT `ghost_nnlo_37`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(1/48*ep^-3 + 1/16*ep^-2 + 7/864*ep^-1)
        //        + cl2*sqrt3*rat(-7/324*ep^-1)
        //        + pi^2*rat(1/192*ep^-1)
        // F1/H = rat(-1/24*ep^-2)
        // F2/H = rat(-1/32*ep^-3 - 13/192*ep^-2 - 23/1152*ep^-1)
        //        + cl2*sqrt3*rat(7/216*ep^-1)
        //        + pi^2*rat(-1/192*ep^-1)
        // F3/H = rat(-1/32*ep^-3 - 3/64*ep^-2 - 49/3456*ep^-1)
        //        + cl2*sqrt3*rat(-7/648*ep^-1)
        //        + pi^2*rat(-1/192*ep^-1)
        // F5/H = rat(1/24*ep^-2)
        // F6/H = rat(1/16*ep^-3 + 1/24*ep^-2)
        //        + pi^2*rat(1/192*ep^-1)
        // F4, F7 = 0
        // sum(Fi)/H = rat(1/48*ep^-3 - 1/96*ep^-2 - 5/192*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[38].renormalization_part(&settings).unwrap();
        // d39: RQFT `ghost_nnlo_38`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-19/384*ep^-3 + 175/768*ep^-2 + 5551/13824*ep^-1)
        //        + cl2*sqrt3*rat(-683/2592*ep^-1)
        //        + pi^2*rat(-19/1536*ep^-1)
        // F1/H = rat(-1/128*ep^-3 - 71/768*ep^-2)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F2/H = rat(5/64*ep^-3 - 77/384*ep^-2 - 149/768*ep^-1)
        //        + cl2*sqrt3*rat(11/144*ep^-1)
        //        + pi^2*rat(5/384*ep^-1)
        // F4/H = rat(9/128*ep^-3 - 73/256*ep^-2 - 1699/13824*ep^-1)
        //        + cl2*sqrt3*rat(485/2592*ep^-1)
        //        + pi^2*rat(3/256*ep^-1)
        // F7/H = rat(1/64*ep^-3 + 17/192*ep^-2)
        //        + pi^2*rat(1/768*ep^-1)
        // F9/H = rat(-5/32*ep^-3 + 7/24*ep^-2)
        //        + pi^2*rat(-5/384*ep^-1)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-19/384*ep^-3 + 23/768*ep^-2 + 65/768*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[39].renormalization_part(&settings).unwrap();
        // d40: RQFT `ghost_nnlo_39`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/384*ep^-3 - 23/768*ep^-2 - 791/13824*ep^-1)
        //        + cl2*sqrt3*rat(103/2592*ep^-1)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F1/H = rat(1/128*ep^-3 + 7/768*ep^-2)
        //        + pi^2*rat(1/1536*ep^-1)
        // F2/H = rat(5/192*ep^-2 + 47/1152*ep^-1)
        //        + cl2*sqrt3*rat(-1/54*ep^-1)
        // F4/H = rat(1/128*ep^-3 + 7/256*ep^-2 + 173/13824*ep^-1)
        //        + cl2*sqrt3*rat(-55/2592*ep^-1)
        //        + pi^2*rat(1/768*ep^-1)
        // F7/H = rat(-1/64*ep^-3 - 1/192*ep^-2)
        //        + pi^2*rat(-1/768*ep^-1)
        // F9/H = rat(-1/48*ep^-2)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-1/384*ep^-3 + 5/768*ep^-2 - 1/256*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[40].renormalization_part(&settings).unwrap();
        // d41: RQFT `ghost_nnlo_40`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/64*ep^-3 - 1/128*ep^-2 - 131/6912*ep^-1)
        //        + cl2*sqrt3*rat(25/1296*ep^-1)
        //        + pi^2*rat(-1/256*ep^-1)
        // F1/H = rat(1/64*ep^-3 - 1/384*ep^-2)
        //        + pi^2*rat(1/768*ep^-1)
        // F2/H = rat(1/64*ep^-3 - 1/384*ep^-2 - 11/2304*ep^-1)
        //        + cl2*sqrt3*rat(1/432*ep^-1)
        //        + pi^2*rat(1/384*ep^-1)
        // F4/H = rat(1/32*ep^-3 + 1/192*ep^-2 + 37/3456*ep^-1)
        //        + cl2*sqrt3*rat(-7/324*ep^-1)
        //        + pi^2*rat(1/192*ep^-1)
        // F7/H = rat(-1/32*ep^-3 + 1/96*ep^-2)
        //        + pi^2*rat(-1/384*ep^-1)
        // F9/H = rat(-1/32*ep^-3 + 1/96*ep^-2)
        //        + pi^2*rat(-1/384*ep^-1)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-1/64*ep^-3 + 5/384*ep^-2 - 5/384*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[41].renormalization_part(&settings).unwrap();
        // d42: RQFT `ghost_nnlo_41`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(1/12*ep^-3 + 9/32*ep^-2 + 35/192*ep^-1)
        //        + cl2*sqrt3*rat(-4/9*ep^-1)
        //        + pi^2*rat(1/48*ep^-1)
        // F1/H = rat(-1/16*ep^-3 - 5/96*ep^-2)
        //        + pi^2*rat(-1/192*ep^-1)
        // F2/H = rat(-3/32*ep^-3 - 13/64*ep^-2 - 23/384*ep^-1)
        //        + cl2*sqrt3*rat(7/72*ep^-1)
        //        + pi^2*rat(-1/64*ep^-1)
        // F3/H = rat(-5/32*ep^-3 - 49/192*ep^-2 - 49/384*ep^-1)
        //        + cl2*sqrt3*rat(25/72*ep^-1)
        //        + pi^2*rat(-5/192*ep^-1)
        // F5/H = rat(1/8*ep^-3)
        //        + pi^2*rat(1/96*ep^-1)
        // F6/H = rat(3/16*ep^-3 + 1/8*ep^-2)
        //        + pi^2*rat(1/64*ep^-1)
        // F4, F7 = 0
        // sum(Fi)/H = rat(1/12*ep^-3 - 5/48*ep^-2 - 1/192*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[42].renormalization_part(&settings).unwrap();
        // d43: RQFT `ghost_nnlo_42`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-79/384*ep^-3 + 323/768*ep^-2 + 659/1536*ep^-1)
        //        + cl2*sqrt3*rat(101/288*ep^-1)
        //        + pi^2*rat(-79/1536*ep^-1)
        // F1/H = rat(19/128*ep^-3 - 307/768*ep^-2)
        //        + pi^2*rat(19/1536*ep^-1)
        // F2/H = rat(15/64*ep^-3 - 77/128*ep^-2 - 149/256*ep^-1)
        //        + cl2*sqrt3*rat(11/48*ep^-1)
        //        + pi^2*rat(5/128*ep^-1)
        // F4/H = rat(49/128*ep^-3 - 427/768*ep^-2 + 173/1536*ep^-1)
        //        + cl2*sqrt3*rat(-167/288*ep^-1)
        //        + pi^2*rat(49/768*ep^-1)
        // F7/H = rat(-19/64*ep^-3 + 37/64*ep^-2)
        //        + pi^2*rat(-19/768*ep^-1)
        // F9/H = rat(-15/32*ep^-3 + 7/8*ep^-2)
        //        + pi^2*rat(-5/128*ep^-1)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-79/384*ep^-3 + 81/256*ep^-2 - 31/768*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[43].renormalization_part(&settings).unwrap();
        // d44: RQFT `ghost_nnlo_43`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/384*ep^-3 - 21/256*ep^-2 - 239/1536*ep^-1)
        //        + cl2*sqrt3*rat(3/32*ep^-1)
        //        + pi^2*rat(-1/1536*ep^-1)
        // F1/H = rat(1/128*ep^-3 + 23/768*ep^-2)
        //        + pi^2*rat(1/1536*ep^-1)
        // F2/H = rat(5/64*ep^-2 + 47/384*ep^-1)
        //        + cl2*sqrt3*rat(-1/18*ep^-1)
        // F4/H = rat(1/128*ep^-3 + 53/768*ep^-2 + 29/1536*ep^-1)
        //        + cl2*sqrt3*rat(-11/288*ep^-1)
        //        + pi^2*rat(1/768*ep^-1)
        // F7/H = rat(-1/64*ep^-3 - 1/64*ep^-2)
        //        + pi^2*rat(-1/768*ep^-1)
        // F9/H = rat(-1/16*ep^-2)
        // F3, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-1/384*ep^-3 + 13/768*ep^-2 - 11/768*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[44].renormalization_part(&settings).unwrap();
        // d45: RQFT `ghost_nnlo_44`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(1/48*ep^-3 + 17/96*ep^-2 + 113/192*ep^-1)
        //        + cl2*sqrt3*rat(-17/36*ep^-1)
        //        + pi^2*rat(1/192*ep^-1)
        // F1/H = rat(-1/16*ep^-3 - 5/96*ep^-2)
        //        + pi^2*rat(-1/192*ep^-1)
        // F2/H = rat(-1/8*ep^-2)
        // F3/H = rat(-1/16*ep^-3 - 23/96*ep^-2 - 33/64*ep^-1)
        //        + cl2*sqrt3*rat(17/36*ep^-1)
        //        + pi^2*rat(-1/96*ep^-1)
        // F5/H = rat(1/8*ep^-3)
        //        + pi^2*rat(1/96*ep^-1)
        // F6/H = rat(1/8*ep^-2)
        // F4, F7 = 0
        // sum(Fi)/H = rat(1/48*ep^-3 - 11/96*ep^-2 + 7/96*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[45].renormalization_part(&settings).unwrap();
        // d46: RQFT `ghost_nnlo_45`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/24*ep^-3 - 89/192*ep^-2 + 277/384*ep^-1)
        //        + cl2*sqrt3*rat(43/36*ep^-1)
        //        + pi^2*rat(-1/96*ep^-1)
        // F1/H = rat(19/128*ep^-3 - 307/768*ep^-2)
        //        + pi^2*rat(19/1536*ep^-1)
        // F3/H = rat(-3/128*ep^-3 - 67/256*ep^-2)
        //        + pi^2*rat(-1/512*ep^-1)
        // F4/H = rat(1/8*ep^-3 + 119/192*ep^-2 - 135/128*ep^-1)
        //        + cl2*sqrt3*rat(-43/36*ep^-1)
        //        + pi^2*rat(1/48*ep^-1)
        // F7/H = rat(-19/64*ep^-3 + 37/64*ep^-2)
        //        + pi^2*rat(-19/768*ep^-1)
        // F9/H = rat(3/64*ep^-3 + 15/64*ep^-2)
        //        + pi^2*rat(1/256*ep^-1)
        // F2, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-1/24*ep^-3 + 59/192*ep^-2 - 1/3*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[46].renormalization_part(&settings).unwrap();
        // d47: RQFT `ghost_nnlo_46`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/96*ep^-3 - 9/64*ep^-1)
        //        + cl2*sqrt3*rat(-1/72*ep^-1)
        //        + pi^2*rat(-1/384*ep^-1)
        // F1/H = rat(1/128*ep^-3 + 23/768*ep^-2)
        //        + pi^2*rat(1/1536*ep^-1)
        // F3/H = rat(3/128*ep^-3 + 3/256*ep^-2)
        //        + pi^2*rat(1/512*ep^-1)
        // F4/H = rat(1/32*ep^-3 - 1/48*ep^-2 + 3/32*ep^-1)
        //        + cl2*sqrt3*rat(1/72*ep^-1)
        //        + pi^2*rat(1/192*ep^-1)
        // F7/H = rat(-1/64*ep^-3 - 1/64*ep^-2)
        //        + pi^2*rat(-1/768*ep^-1)
        // F9/H = rat(-3/64*ep^-3 + 1/64*ep^-2)
        //        + pi^2*rat(-1/256*ep^-1)
        // F2, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-1/96*ep^-3 + 1/48*ep^-2 - 3/64*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[47].renormalization_part(&settings).unwrap();
        // d48: RQFT `ghost_nnlo_47`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/96*ep^-3 - 1/64*ep^-2 + 1/384*ep^-1)
        //        + cl2*sqrt3*rat(1/72*ep^-1)
        //        + pi^2*rat(-1/384*ep^-1)
        // F1/H = rat(1/64*ep^-3 - 1/384*ep^-2)
        //        + pi^2*rat(1/768*ep^-1)
        // F3/H = rat(1/64*ep^-3 - 1/384*ep^-2)
        //        + pi^2*rat(1/768*ep^-1)
        // F4/H = rat(1/32*ep^-3 + 1/64*ep^-2 - 5/384*ep^-1)
        //        + cl2*sqrt3*rat(-1/72*ep^-1)
        //        + pi^2*rat(1/192*ep^-1)
        // F7/H = rat(-1/32*ep^-3 + 1/96*ep^-2)
        //        + pi^2*rat(-1/384*ep^-1)
        // F9/H = rat(-1/32*ep^-3 + 1/96*ep^-2)
        //        + pi^2*rat(-1/384*ep^-1)
        // F2, F5..F6, F8, F10..F11 = 0
        // sum(Fi)/H = rat(-1/96*ep^-3 + 1/64*ep^-2 - 1/96*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[48].renormalization_part(&settings).unwrap();
        // d49: RQFT `ghost_nnlo_48`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/32*ep^-3 - 5/192*ep^-2 - 5/384*ep^-1)
        //        + cl2*sqrt3*rat(-1/72*ep^-1)
        //        + z3*rat(1/8*ep^-1)
        //        + pi^2*rat(-1/128*ep^-1)
        // F1/H = rat(3/64*ep^-3 - 1/128*ep^-2 - 11/768*ep^-1)
        //        + cl2*sqrt3*rat(1/144*ep^-1)
        //        + pi^2*rat(1/128*ep^-1)
        // F5/H = rat(3/32*ep^-3 + 3/64*ep^-2 - 5/192*ep^-1)
        //        + pi^2*rat(1/128*ep^-1)
        // F8/H = rat(3/64*ep^-3 - 1/128*ep^-2 - 11/768*ep^-1)
        //        + cl2*sqrt3*rat(1/144*ep^-1)
        //        + pi^2*rat(1/128*ep^-1)
        // F12/H = rat(-3/32*ep^-3 + 1/32*ep^-2)
        //         + pi^2*rat(-1/128*ep^-1)
        // F15/H = rat(-3/32*ep^-3 + 1/32*ep^-2)
        //         + pi^2*rat(-1/128*ep^-1)
        // F2..F4, F6..F7, F9..F11, F13..F14, F16..F17 = 0
        // sum(Fi)/H = rat(-1/32*ep^-3 + 13/192*ep^-2 - 13/192*ep^-1)
        //             + z3*rat(1/8*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[49].renormalization_part(&settings).unwrap();
        // d50: RQFT `ghost_nnlo_49`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/96*ep^-3 - 1/192*ep^-2 + 5/1152*ep^-1)
        //        + cl2*sqrt3*rat(-1/216*ep^-1)
        //        + pi^2*rat(-1/384*ep^-1)
        // F1/H = rat(1/64*ep^-3 - 1/384*ep^-2 - 11/2304*ep^-1)
        //        + cl2*sqrt3*rat(1/432*ep^-1)
        //        + pi^2*rat(1/384*ep^-1)
        // F5/H = rat(1/32*ep^-3 + 1/192*ep^-2 - 1/192*ep^-1)
        //        + pi^2*rat(1/384*ep^-1)
        // F8/H = rat(1/64*ep^-3 - 1/384*ep^-2 - 11/2304*ep^-1)
        //        + cl2*sqrt3*rat(1/432*ep^-1)
        //        + pi^2*rat(1/384*ep^-1)
        // F12/H = rat(-1/32*ep^-3 + 1/96*ep^-2)
        //         + pi^2*rat(-1/384*ep^-1)
        // F15/H = rat(-1/32*ep^-3 + 1/96*ep^-2)
        //         + pi^2*rat(-1/384*ep^-1)
        // F2..F4, F6..F7, F9..F11, F13..F14, F16..F17 = 0
        // sum(Fi)/H = rat(-1/96*ep^-3 + 1/64*ep^-2 - 1/96*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[50].renormalization_part(&settings).unwrap();
        // d51: RQFT `ghost_nnlo_50`.
        // H = p1.p1*gs^6*ca*nf
        // F0/H = cf*rat(1/6*ep^-3 + 29/18*ep^-2 + 59/18*ep^-1)
        //        + cf*cl2*sqrt3*rat(-14/27*ep^-1)
        //        + ca*rat(-1/12*ep^-3 - 29/36*ep^-2 - 59/36*ep^-1)
        //        + ca*cl2*sqrt3*rat(7/27*ep^-1)
        //        + z3*cf*rat(-2*ep^-1)
        //        + z3*ca*rat(ep^-1)
        //        + pi^2*cf*rat(1/24*ep^-1)
        //        + pi^2*ca*rat(-1/48*ep^-1)
        // F1/H = cf*rat(-1/4*ep^-3 - 13/24*ep^-2 - 23/144*ep^-1)
        //        + cf*cl2*sqrt3*rat(7/27*ep^-1)
        //        + ca*rat(1/8*ep^-3 + 13/48*ep^-2 + 23/288*ep^-1)
        //        + ca*cl2*sqrt3*rat(-7/54*ep^-1)
        //        + pi^2*cf*rat(-1/24*ep^-1)
        //        + pi^2*ca*rat(1/48*ep^-1)
        // F2/H = cf*rat(-4/3*ep^-2 - 4/9*ep^-1)
        //        + ca*rat(2/3*ep^-2 + 2/9*ep^-1)
        // F3/H = cf*rat(-1/2*ep^-3 - 5/3*ep^-2 - 8/9*ep^-1)
        //        + ca*rat(1/4*ep^-3 + 5/6*ep^-2 + 4/9*ep^-1)
        //        + pi^2*cf*rat(-1/24*ep^-1)
        //        + pi^2*ca*rat(1/48*ep^-1)
        // F4/H = cf*rat(-1/4*ep^-3 - 13/24*ep^-2 - 23/144*ep^-1)
        //        + cf*cl2*sqrt3*rat(7/27*ep^-1)
        //        + ca*rat(1/8*ep^-3 + 13/48*ep^-2 + 23/288*ep^-1)
        //        + ca*cl2*sqrt3*rat(-7/54*ep^-1)
        //        + pi^2*cf*rat(-1/24*ep^-1)
        //        + pi^2*ca*rat(1/48*ep^-1)
        // F6/H = cf*rat(1/2*ep^-3 + 1/6*ep^-2 + 1/9*ep^-1)
        //        + ca*rat(-1/4*ep^-3 - 1/12*ep^-2 - 1/18*ep^-1)
        //        + pi^2*cf*rat(1/24*ep^-1)
        //        + pi^2*ca*rat(-1/48*ep^-1)
        // F7/H = cf*rat(4/3*ep^-2 + 4/9*ep^-1)
        //        + ca*rat(-2/3*ep^-2 - 2/9*ep^-1)
        // F9/H = cf*rat(1/2*ep^-3 + 1/6*ep^-2 + 1/9*ep^-1)
        //        + ca*rat(-1/4*ep^-3 - 1/12*ep^-2 - 1/18*ep^-1)
        //        + pi^2*cf*rat(1/24*ep^-1)
        //        + pi^2*ca*rat(-1/48*ep^-1)
        // F5, F8 = 0
        // sum(Fi)/H = cf*rat(1/6*ep^-3 - 29/36*ep^-2 + 55/24*ep^-1)
        //             + ca*rat(-1/12*ep^-3 + 29/72*ep^-2 - 55/48*ep^-1)
        //             + z3*cf*rat(-2*ep^-1)
        //             + z3*ca*rat(ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[51].renormalization_part(&settings).unwrap();
        // d52: RQFT `ghost_nnlo_51`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(-1/72*ep^-3 + 155/144*ep^-2 + 751/288*ep^-1)
        //        + cl2*sqrt3*rat(-43/54*ep^-1)
        //        + pi^2*rat(-1/288*ep^-1)
        // F1/H = rat(5/12*ep^-3 - 77/72*ep^-2 - 149/144*ep^-1)
        //        + cl2*sqrt3*rat(11/27*ep^-1)
        //        + pi^2*rat(5/72*ep^-1)
        // F4/H = rat(1/24*ep^-3 - 151/144*ep^-2 - 95/72*ep^-1)
        //        + pi^2*rat(1/288*ep^-1)
        // F5/H = rat(-3/8*ep^-3 - 13/16*ep^-2 - 23/96*ep^-1)
        //        + cl2*sqrt3*rat(7/18*ep^-1)
        //        + pi^2*rat(-1/16*ep^-1)
        // F8/H = rat(-5/6*ep^-3 + 4/3*ep^-2 + 13/27*ep^-1)
        //        + pi^2*rat(-5/72*ep^-1)
        // F9/H = rat(3/4*ep^-3 + 1/4*ep^-2 + 1/6*ep^-1)
        //        + pi^2*rat(1/16*ep^-1)
        // F2..F3, F6..F7 = 0
        // sum(Fi)/H = rat(-1/72*ep^-3 - 13/48*ep^-2 + 143/216*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[52].renormalization_part(&settings).unwrap();
        // d53: RQFT `ghost_nnlo_52`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(-1/72*ep^-3 + 155/144*ep^-2 + 751/288*ep^-1)
        //        + cl2*sqrt3*rat(-43/54*ep^-1)
        //        + pi^2*rat(-1/288*ep^-1)
        // F1/H = rat(-3/8*ep^-3 - 13/16*ep^-2 - 23/96*ep^-1)
        //        + cl2*sqrt3*rat(7/18*ep^-1)
        //        + pi^2*rat(-1/16*ep^-1)
        // F2/H = rat(1/24*ep^-3 - 151/144*ep^-2 - 95/72*ep^-1)
        //        + pi^2*rat(1/288*ep^-1)
        // F5/H = rat(5/12*ep^-3 - 77/72*ep^-2 - 149/144*ep^-1)
        //        + cl2*sqrt3*rat(11/27*ep^-1)
        //        + pi^2*rat(5/72*ep^-1)
        // F6/H = rat(3/4*ep^-3 + 1/4*ep^-2 + 1/6*ep^-1)
        //        + pi^2*rat(1/16*ep^-1)
        // F7/H = rat(-5/6*ep^-3 + 4/3*ep^-2 + 13/27*ep^-1)
        //        + pi^2*rat(-5/72*ep^-1)
        // F3..F4, F8..F9 = 0
        // sum(Fi)/H = rat(-1/72*ep^-3 - 13/48*ep^-2 + 143/216*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[53].renormalization_part(&settings).unwrap();
        // d54: RQFT `ghost_nnlo_53`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-65/96*ep^-3 + 3329/1152*ep^-2 + 23713/2304*ep^-1)
        //        + cl2*sqrt3*rat(-143/72*ep^-1)
        //        + z3*rat(1/4*ep^-1)
        //        + pi^2*rat(-65/384*ep^-1)
        // F1/H = rat(65/64*ep^-3 - 1001/384*ep^-2 - 1937/768*ep^-1)
        //        + cl2*sqrt3*rat(143/144*ep^-1)
        //        + pi^2*rat(65/384*ep^-1)
        // F3/H = rat(-9/8*ep^-2 - 8/9*ep^-1)
        // F5/H = rat(65/32*ep^-3 - 463/128*ep^-2 - 3685/576*ep^-1)
        //        + pi^2*rat(65/384*ep^-1)
        // F8/H = rat(65/64*ep^-3 - 1001/384*ep^-2 - 1937/768*ep^-1)
        //        + cl2*sqrt3*rat(143/144*ep^-1)
        //        + pi^2*rat(65/384*ep^-1)
        // F12/H = rat(-65/32*ep^-3 + 13/4*ep^-2 + 169/144*ep^-1)
        //         + pi^2*rat(-65/384*ep^-1)
        // F13/H = rat(9/8*ep^-2 + 8/9*ep^-1)
        // F15/H = rat(-65/32*ep^-3 + 13/4*ep^-2 + 169/144*ep^-1)
        //         + pi^2*rat(-65/384*ep^-1)
        // F2, F4, F6..F7, F9..F11, F14, F16..F17 = 0
        // sum(Fi)/H = rat(-65/96*ep^-3 + 161/288*ep^-2 + 2759/2304*ep^-1)
        //             + z3*rat(1/4*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[54].renormalization_part(&settings).unwrap();
        // d55: RQFT `ghost_nnlo_54`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(5/576*ep^-3 - 13/192*ep^-2 - 319/1152*ep^-1)
        //        + cl2*sqrt3*rat(35/432*ep^-1)
        //        + pi^2*rat(5/2304*ep^-1)
        // F1/H = rat(5/64*ep^-2 + 47/384*ep^-1)
        //        + cl2*sqrt3*rat(-1/18*ep^-1)
        // F5/H = rat(-5/192*ep^-3 + 17/288*ep^-2 + 23/144*ep^-1)
        //        + pi^2*rat(-5/2304*ep^-1)
        // F8/H = rat(-5/192*ep^-3 + 77/1152*ep^-2 + 149/2304*ep^-1)
        //        + cl2*sqrt3*rat(-11/432*ep^-1)
        //        + pi^2*rat(-5/1152*ep^-1)
        // F12/H = rat(-1/16*ep^-2 - 1/16*ep^-1)
        // F15/H = rat(5/96*ep^-3 - 1/12*ep^-2 - 13/432*ep^-1)
        //         + pi^2*rat(5/1152*ep^-1)
        // F2..F4, F6..F7, F9..F11, F13..F14, F16..F17 = 0
        // sum(Fi)/H = rat(5/576*ep^-3 - 11/1152*ep^-2 - 157/6912*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[55].renormalization_part(&settings).unwrap();
        // d56: RQFT `ghost_nnlo_55`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(5/576*ep^-3 - 13/192*ep^-2 - 319/1152*ep^-1)
        //        + cl2*sqrt3*rat(35/432*ep^-1)
        //        + pi^2*rat(5/2304*ep^-1)
        // F1/H = rat(-5/192*ep^-3 + 77/1152*ep^-2 + 149/2304*ep^-1)
        //        + cl2*sqrt3*rat(-11/432*ep^-1)
        //        + pi^2*rat(-5/1152*ep^-1)
        // F5/H = rat(-5/192*ep^-3 + 17/288*ep^-2 + 23/144*ep^-1)
        //        + pi^2*rat(-5/2304*ep^-1)
        // F8/H = rat(5/64*ep^-2 + 47/384*ep^-1)
        //        + cl2*sqrt3*rat(-1/18*ep^-1)
        // F12/H = rat(5/96*ep^-3 - 1/12*ep^-2 - 13/432*ep^-1)
        //         + pi^2*rat(5/1152*ep^-1)
        // F15/H = rat(-1/16*ep^-2 - 1/16*ep^-1)
        // F2..F4, F6..F7, F9..F11, F13..F14, F16..F17 = 0
        // sum(Fi)/H = rat(5/576*ep^-3 - 11/1152*ep^-2 - 157/6912*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[56].renormalization_part(&settings).unwrap();
        // d57: RQFT `ghost_nnlo_56`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-59/1152*ep^-2 - 161/768*ep^-1)
        //        + cl2*sqrt3*rat(1/27*ep^-1)
        // F1/H = rat(5/192*ep^-2 + 47/1152*ep^-1)
        //        + cl2*sqrt3*rat(-1/54*ep^-1)
        // F3/H = rat(1/24*ep^-2 + 1/36*ep^-1)
        // F5/H = rat(19/384*ep^-2 + 25/192*ep^-1)
        // F8/H = rat(5/192*ep^-2 + 47/1152*ep^-1)
        //        + cl2*sqrt3*rat(-1/54*ep^-1)
        // F12/H = rat(-1/48*ep^-2 - 1/48*ep^-1)
        // F13/H = rat(-1/24*ep^-2 - 1/36*ep^-1)
        // F15/H = rat(-1/48*ep^-2 - 1/48*ep^-1)
        // F2, F4, F6..F7, F9..F11, F14, F16..F17 = 0
        // sum(Fi)/H = rat(5/576*ep^-2 - 91/2304*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[57].renormalization_part(&settings).unwrap();
        // d58: RQFT `ghost_nnlo_57`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/24*ep^-3 + 1/32*ep^-2 + 29/216*ep^-1)
        //        + cl2*sqrt3*rat(-13/162*ep^-1)
        //        + pi^2*rat(-1/96*ep^-1)
        // F5/H = rat(1/16*ep^-3 - 1/24*ep^-2 - 31/864*ep^-1)
        //        + cl2*sqrt3*rat(13/324*ep^-1)
        //        + pi^2*rat(1/96*ep^-1)
        // F6/H = rat(1/16*ep^-3 - 1/24*ep^-2 - 31/864*ep^-1)
        //        + cl2*sqrt3*rat(13/324*ep^-1)
        //        + pi^2*rat(1/96*ep^-1)
        // F11/H = rat(-1/8*ep^-3 + 7/96*ep^-2 - 1/96*ep^-1)
        //         + pi^2*rat(-1/96*ep^-1)
        // F1..F4, F7..F10 = 0
        // sum(Fi)/H = rat(-1/24*ep^-3 + 1/48*ep^-2 + 5/96*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[58].renormalization_part(&settings).unwrap();
        // d59: RQFT `ghost_nnlo_58`.
        // H = p1.p1*gs^6*ca*nf^2
        // F0/H = rat(-1/9*ep^-3 - 55/54*ep^-2 - 161/162*ep^-1)
        //        + cl2*sqrt3*rat(92/81*ep^-1)
        //        + pi^2*rat(-1/36*ep^-1)
        // F1/H = rat(1/6*ep^-3 + 35/36*ep^-2 + 133/216*ep^-1)
        //        + cl2*sqrt3*rat(-46/81*ep^-1)
        //        + pi^2*rat(1/36*ep^-1)
        // F2/H = rat(1/6*ep^-3 + 35/36*ep^-2 + 133/216*ep^-1)
        //        + cl2*sqrt3*rat(-46/81*ep^-1)
        //        + pi^2*rat(1/36*ep^-1)
        // F3/H = rat(-1/3*ep^-3 - 5/6*ep^-2 - 7/54*ep^-1)
        //        + pi^2*rat(-1/36*ep^-1)
        // sum(Fi)/H = rat(-1/9*ep^-3 + 5/54*ep^-2 + 35/324*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[59].renormalization_part(&settings).unwrap();
        // d60: RQFT `ghost_nnlo_59`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(19/72*ep^-3 - 569/432*ep^-2 - 2489/648*ep^-1)
        //        + cl2*sqrt3*rat(203/162*ep^-1)
        //        + pi^2*rat(19/288*ep^-1)
        // F3/H = rat(-19/48*ep^-3 + 397/288*ep^-2 + 1429/576*ep^-1)
        //        + cl2*sqrt3*rat(-37/108*ep^-1)
        //        + pi^2*rat(-19/288*ep^-1)
        // F4/H = rat(-19/48*ep^-3 + 439/288*ep^-2 + 3617/1728*ep^-1)
        //        + cl2*sqrt3*rat(-295/324*ep^-1)
        //        + pi^2*rat(-19/288*ep^-1)
        // F7/H = rat(19/24*ep^-3 - 89/48*ep^-2 - 431/432*ep^-1)
        //        + pi^2*rat(19/288*ep^-1)
        // F1..F2, F5..F6 = 0
        // sum(Fi)/H = rat(19/72*ep^-3 - 29/108*ep^-2 - 343/1296*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[60].renormalization_part(&settings).unwrap();
        // d61: RQFT `ghost_nnlo_60`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(1/72*ep^-3 + 97/432*ep^-2 + 179/324*ep^-1)
        //        + cl2*sqrt3*rat(-13/54*ep^-1)
        //        + pi^2*rat(1/288*ep^-1)
        // F3/H = rat(-1/48*ep^-3 - 65/288*ep^-2 - 803/1728*ep^-1)
        //        + cl2*sqrt3*rat(43/324*ep^-1)
        //        + pi^2*rat(-1/288*ep^-1)
        // F4/H = rat(-1/48*ep^-3 - 59/288*ep^-2 - 301/1728*ep^-1)
        //        + cl2*sqrt3*rat(35/324*ep^-1)
        //        + pi^2*rat(-1/288*ep^-1)
        // F7/H = rat(1/24*ep^-3 + 3/16*ep^-2 + 31/432*ep^-1)
        //        + pi^2*rat(1/288*ep^-1)
        // F1..F2, F5..F6 = 0
        // sum(Fi)/H = rat(1/72*ep^-3 - 1/54*ep^-2 - 19/1296*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[61].renormalization_part(&settings).unwrap();
        // d62: RQFT `ghost_nnlo_61`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(19/72*ep^-3 - 569/432*ep^-2 - 2489/648*ep^-1)
        //        + cl2*sqrt3*rat(203/162*ep^-1)
        //        + pi^2*rat(19/288*ep^-1)
        // F3/H = rat(-19/48*ep^-3 + 439/288*ep^-2 + 3617/1728*ep^-1)
        //        + cl2*sqrt3*rat(-295/324*ep^-1)
        //        + pi^2*rat(-19/288*ep^-1)
        // F4/H = rat(-19/48*ep^-3 + 397/288*ep^-2 + 1429/576*ep^-1)
        //        + cl2*sqrt3*rat(-37/108*ep^-1)
        //        + pi^2*rat(-19/288*ep^-1)
        // F7/H = rat(19/24*ep^-3 - 89/48*ep^-2 - 431/432*ep^-1)
        //        + pi^2*rat(19/288*ep^-1)
        // F1..F2, F5..F6 = 0
        // sum(Fi)/H = rat(19/72*ep^-3 - 29/108*ep^-2 - 343/1296*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[62].renormalization_part(&settings).unwrap();
        // d63: RQFT `ghost_nnlo_62`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-179/288*ep^-3 + 5153/3456*ep^-2 - 17635/5184*ep^-1)
        //        + cl2*sqrt3*rat(-575/216*ep^-1)
        //        + pi^2*rat(-179/1152*ep^-1)
        // F5/H = rat(179/192*ep^-3 - 1075/576*ep^-2 + 5399/1152*ep^-1)
        //        + cl2*sqrt3*rat(575/432*ep^-1)
        //        + pi^2*rat(179/1152*ep^-1)
        // F6/H = rat(179/192*ep^-3 - 1075/576*ep^-2 + 5399/1152*ep^-1)
        //        + cl2*sqrt3*rat(575/432*ep^-1)
        //        + pi^2*rat(179/1152*ep^-1)
        // F11/H = rat(-179/96*ep^-3 + 383/128*ep^-2 - 18715/3456*ep^-1)
        //         + pi^2*rat(-179/1152*ep^-1)
        // F1..F4, F7..F10 = 0
        // sum(Fi)/H = rat(-179/288*ep^-3 + 1297/1728*ep^-2 + 5767/10368*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[63].renormalization_part(&settings).unwrap();
        // d64: RQFT `ghost_nnlo_63`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-11/288*ep^-3 + 455/3456*ep^-2 + 2629/2592*ep^-1)
        //        + cl2*sqrt3*rat(-77/648*ep^-1)
        //        + pi^2*rat(-11/1152*ep^-1)
        // F5/H = rat(11/192*ep^-3 - 41/288*ep^-2 - 1277/1728*ep^-1)
        //        + cl2*sqrt3*rat(121/1296*ep^-1)
        //        + pi^2*rat(11/1152*ep^-1)
        // F6/H = rat(11/192*ep^-3 - 25/144*ep^-2 - 359/576*ep^-1)
        //        + cl2*sqrt3*rat(11/432*ep^-1)
        //        + pi^2*rat(11/1152*ep^-1)
        // F11/H = rat(-11/96*ep^-3 + 91/384*ep^-2 + 1307/3456*ep^-1)
        //         + pi^2*rat(-11/1152*ep^-1)
        // F1..F4, F7..F10 = 0
        // sum(Fi)/H = rat(-11/288*ep^-3 + 91/1728*ep^-2 + 313/10368*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[64].renormalization_part(&settings).unwrap();
        // d65: RQFT `ghost_nnlo_64`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(1/72*ep^-3 + 97/432*ep^-2 + 179/324*ep^-1)
        //        + cl2*sqrt3*rat(-13/54*ep^-1)
        //        + pi^2*rat(1/288*ep^-1)
        // F3/H = rat(-1/48*ep^-3 - 59/288*ep^-2 - 301/1728*ep^-1)
        //        + cl2*sqrt3*rat(35/324*ep^-1)
        //        + pi^2*rat(-1/288*ep^-1)
        // F4/H = rat(-1/48*ep^-3 - 65/288*ep^-2 - 803/1728*ep^-1)
        //        + cl2*sqrt3*rat(43/324*ep^-1)
        //        + pi^2*rat(-1/288*ep^-1)
        // F7/H = rat(1/24*ep^-3 + 3/16*ep^-2 + 31/432*ep^-1)
        //        + pi^2*rat(1/288*ep^-1)
        // F1..F2, F5..F6 = 0
        // sum(Fi)/H = rat(1/72*ep^-3 - 1/54*ep^-2 - 19/1296*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[65].renormalization_part(&settings).unwrap();
        // d66: RQFT `ghost_nnlo_65`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-11/288*ep^-3 + 455/3456*ep^-2 + 2629/2592*ep^-1)
        //        + cl2*sqrt3*rat(-77/648*ep^-1)
        //        + pi^2*rat(-11/1152*ep^-1)
        // F5/H = rat(11/192*ep^-3 - 25/144*ep^-2 - 359/576*ep^-1)
        //        + cl2*sqrt3*rat(11/432*ep^-1)
        //        + pi^2*rat(11/1152*ep^-1)
        // F6/H = rat(11/192*ep^-3 - 41/288*ep^-2 - 1277/1728*ep^-1)
        //        + cl2*sqrt3*rat(121/1296*ep^-1)
        //        + pi^2*rat(11/1152*ep^-1)
        // F11/H = rat(-11/96*ep^-3 + 91/384*ep^-2 + 1307/3456*ep^-1)
        //         + pi^2*rat(-11/1152*ep^-1)
        // F1..F4, F7..F10 = 0
        // sum(Fi)/H = rat(-11/288*ep^-3 + 91/1728*ep^-2 + 313/10368*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[66].renormalization_part(&settings).unwrap();
        // d67: RQFT `ghost_nnlo_66`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(1/288*ep^-3 - 175/3456*ep^-2 - 1081/5184*ep^-1)
        //        + cl2*sqrt3*rat(47/648*ep^-1)
        //        + pi^2*rat(1/1152*ep^-1)
        // F5/H = rat(-1/192*ep^-3 + 29/576*ep^-2 + 395/3456*ep^-1)
        //        + cl2*sqrt3*rat(-47/1296*ep^-1)
        //        + pi^2*rat(-1/1152*ep^-1)
        // F6/H = rat(-1/192*ep^-3 + 29/576*ep^-2 + 395/3456*ep^-1)
        //        + cl2*sqrt3*rat(-47/1296*ep^-1)
        //        + pi^2*rat(-1/1152*ep^-1)
        // F11/H = rat(1/96*ep^-3 - 19/384*ep^-2 - 91/3456*ep^-1)
        //         + pi^2*rat(1/1152*ep^-1)
        // F1..F4, F7..F10 = 0
        // sum(Fi)/H = rat(1/288*ep^-3 + 1/1728*ep^-2 - 65/10368*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[67].renormalization_part(&settings).unwrap();
        // d68: RQFT `ghost_nnlo_67`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/48*ep^-3 - 1/96*ep^-2 - 1/1728*ep^-1)
        //        + cl2*sqrt3*rat(5/324*ep^-1)
        //        + pi^2*rat(-1/192*ep^-1)
        // F3/H = rat(1/16*ep^-3 - 1/32*ep^-2 + 1/288*ep^-1)
        //        + pi^2*rat(1/192*ep^-1)
        // F5/H = rat(1/16*ep^-3 - 1/96*ep^-2 - 17/1728*ep^-1)
        //        + cl2*sqrt3*rat(-5/324*ep^-1)
        //        + pi^2*rat(1/96*ep^-1)
        // F13/H = rat(-1/8*ep^-3 + 1/12*ep^-2 - 1/72*ep^-1)
        //         + pi^2*rat(-1/96*ep^-1)
        // F1..F2, F4, F6..F12, F14..F17 = 0
        // sum(Fi)/H = rat(-1/48*ep^-3 + 1/32*ep^-2 - 1/48*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[68].renormalization_part(&settings).unwrap();
        // d69: RQFT `ghost_nnlo_68`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(1/24*ep^-3 - 1/48*ep^-2 - 233/864*ep^-1)
        //        + cl2*sqrt3*rat(67/162*ep^-1)
        //        + pi^2*rat(1/96*ep^-1)
        // F1/H = rat(-1/8*ep^-3 - 11/48*ep^-2 + 13/144*ep^-1)
        //        + pi^2*rat(-1/96*ep^-1)
        // F3/H = rat(-1/8*ep^-3 + 5/48*ep^-2 + 293/864*ep^-1)
        //        + cl2*sqrt3*rat(-67/162*ep^-1)
        //        + pi^2*rat(-1/48*ep^-1)
        // F5/H = rat(1/4*ep^-3 + 1/12*ep^-2 - 1/18*ep^-1)
        //        + pi^2*rat(1/48*ep^-1)
        // F2, F4, F6..F7 = 0
        // sum(Fi)/H = rat(1/24*ep^-3 - 1/16*ep^-2 + 5/48*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[69].renormalization_part(&settings).unwrap();
        // d70: RQFT `ghost_nnlo_69`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-5/48*ep^-3 + 1/32*ep^-2 - 1369/1728*ep^-1)
        //        + cl2*sqrt3*rat(257/324*ep^-1)
        //        + pi^2*rat(-5/192*ep^-1)
        // F3/H = rat(5/16*ep^-3 - 29/32*ep^-2 + 77/288*ep^-1)
        //        + pi^2*rat(5/192*ep^-1)
        // F5/H = rat(5/16*ep^-3 - 9/32*ep^-2 + 1057/1728*ep^-1)
        //        + cl2*sqrt3*rat(-257/324*ep^-1)
        //        + pi^2*rat(5/96*ep^-1)
        // F13/H = rat(-5/8*ep^-3 + 11/8*ep^-2 - 7/18*ep^-1)
        //         + pi^2*rat(-5/96*ep^-1)
        // F1..F2, F4, F6..F12, F14..F17 = 0
        // sum(Fi)/H = rat(-5/48*ep^-3 + 7/32*ep^-2 - 29/96*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[70].renormalization_part(&settings).unwrap();
        // d71: RQFT `ghost_nnlo_70`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-1/48*ep^-2 + 25/864*ep^-1)
        //        + cl2*sqrt3*rat(-4/81*ep^-1)
        // F3/H = rat(5/48*ep^-2 - 5/144*ep^-1)
        // F5/H = rat(1/48*ep^-2 - 37/864*ep^-1)
        //        + cl2*sqrt3*rat(4/81*ep^-1)
        // F13/H = rat(-1/12*ep^-2 + 1/36*ep^-1)
        // F1..F2, F4, F6..F12, F14..F17 = 0
        // sum(Fi)/H = rat(1/48*ep^-2 - 1/48*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[71].renormalization_part(&settings).unwrap();
        // d72: RQFT `ghost_nnlo_71`.
        // H = p1.p1*gs^6*ca*cf*nf
        // F0/H = rat(-1/12*ep^-3 - 61/72*ep^-2 - 431/432*ep^-1)
        //        + cl2*sqrt3*rat(49/81*ep^-1)
        //        + pi^2*rat(-1/48*ep^-1)
        // F1/H = rat(2/3*ep^-2 + 2/9*ep^-1)
        // F2/H = rat(1/4*ep^-3 + 3/8*ep^-2 + 11/18*ep^-1)
        //        + pi^2*rat(1/48*ep^-1)
        // F3/H = rat(1/4*ep^-3 + 7/8*ep^-2 + 53/432*ep^-1)
        //        + cl2*sqrt3*rat(-49/81*ep^-1)
        //        + pi^2*rat(1/24*ep^-1)
        // F5/H = rat(-2/3*ep^-2 - 2/9*ep^-1)
        // F7/H = rat(-1/2*ep^-3 - 1/4*ep^-2 + 1/18*ep^-1)
        //        + pi^2*rat(-1/24*ep^-1)
        // F4, F6 = 0
        // sum(Fi)/H = rat(-1/12*ep^-3 + 11/72*ep^-2 - 5/24*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[72].renormalization_part(&settings).unwrap();
        // d73: RQFT `ghost_nnlo_72`.
        // H = p1.p1*gs^6*ca*cf*nf
        // F0/H = rat(-1/12*ep^-3 - 61/72*ep^-2 - 431/432*ep^-1)
        //        + cl2*sqrt3*rat(49/81*ep^-1)
        //        + pi^2*rat(-1/48*ep^-1)
        // F1/H = rat(2/3*ep^-2 + 2/9*ep^-1)
        // F2/H = rat(1/4*ep^-3 + 3/8*ep^-2 + 11/18*ep^-1)
        //        + pi^2*rat(1/48*ep^-1)
        // F3/H = rat(1/4*ep^-3 + 7/8*ep^-2 + 53/432*ep^-1)
        //        + cl2*sqrt3*rat(-49/81*ep^-1)
        //        + pi^2*rat(1/24*ep^-1)
        // F5/H = rat(-2/3*ep^-2 - 2/9*ep^-1)
        // F7/H = rat(-1/2*ep^-3 - 1/4*ep^-2 + 1/18*ep^-1)
        //        + pi^2*rat(-1/24*ep^-1)
        // F4, F6 = 0
        // sum(Fi)/H = rat(-1/12*ep^-3 + 11/72*ep^-2 - 5/24*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[73].renormalization_part(&settings).unwrap();
        // d74: RQFT `ghost_nnlo_73`.
        // H = p1.p1*gs^6*ca^2*nf
        // F0/H = rat(23/72*ep^-3 + 25/144*ep^-2 - 373/96*ep^-1)
        //        + cl2*sqrt3*rat(85/54*ep^-1)
        //        + pi^2*rat(23/288*ep^-1)
        // F1/H = rat(-23/24*ep^-3 - 253/144*ep^-2 + 47/9*ep^-1)
        //        + pi^2*rat(-23/288*ep^-1)
        // F3/H = rat(-23/24*ep^-3 + 23/144*ep^-2 + 221/96*ep^-1)
        //        + cl2*sqrt3*rat(-85/54*ep^-1)
        //        + pi^2*rat(-23/144*ep^-1)
        // F5/H = rat(23/12*ep^-3 + 11/12*ep^-2 - 169/54*ep^-1)
        //        + pi^2*rat(23/144*ep^-1)
        // F2, F4, F6..F7 = 0
        // sum(Fi)/H = rat(23/72*ep^-3 - 73/144*ep^-2 + 55/108*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[74].renormalization_part(&settings).unwrap();
        // d75: RQFT `ghost_nnlo_74`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-223/288*ep^-3 + 5615/576*ep^-2 + 20179/1152*ep^-1)
        //        + cl2*sqrt3*rat(-401/216*ep^-1)
        //        + pi^2*rat(-223/1152*ep^-1)
        // F1/H = rat(-7/2*ep^-2 - 79/36*ep^-1)
        // F3/H = rat(223/96*ep^-3 - 9799/576*ep^-2 - 1477/288*ep^-1)
        //        + pi^2*rat(223/1152*ep^-1)
        // F5/H = rat(223/96*ep^-3 - 8215/576*ep^-2 - 607/128*ep^-1)
        //        + cl2*sqrt3*rat(401/216*ep^-1)
        //        + pi^2*rat(223/576*ep^-1)
        // F8/H = rat(-7/2*ep^-2 - 79/36*ep^-1)
        // F9/H = rat(7/2*ep^-2 + 79/36*ep^-1)
        // F13/H = rat(-223/48*ep^-3 + 559/24*ep^-2 - 260/27*ep^-1)
        //         + pi^2*rat(-223/576*ep^-1)
        // F14/H = rat(7/2*ep^-2 + 79/36*ep^-1)
        // F2, F4, F6..F7, F10..F12, F15..F17 = 0
        // sum(Fi)/H = rat(-223/288*ep^-3 + 113/64*ep^-2 - 857/432*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[75].renormalization_part(&settings).unwrap();
        // d76: RQFT `ghost_nnlo_75`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-7/288*ep^-3 - 145/576*ep^-2 + 283/1152*ep^-1)
        //        + cl2*sqrt3*rat(-17/216*ep^-1)
        //        + pi^2*rat(-7/1152*ep^-1)
        // F3/H = rat(7/96*ep^-3 + 497/576*ep^-2 - 319/288*ep^-1)
        //        + pi^2*rat(7/1152*ep^-1)
        // F5/H = rat(7/96*ep^-3 + 209/576*ep^-2 - 101/384*ep^-1)
        //        + cl2*sqrt3*rat(17/216*ep^-1)
        //        + pi^2*rat(7/576*ep^-1)
        // F13/H = rat(-7/48*ep^-3 - 5/6*ep^-2 + 103/108*ep^-1)
        //         + pi^2*rat(-7/576*ep^-1)
        // F1..F2, F4, F6..F12, F14..F17 = 0
        // sum(Fi)/H = rat(-7/288*ep^-3 + 9/64*ep^-2 - 37/216*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[76].renormalization_part(&settings).unwrap();
        // d77: RQFT `ghost_nnlo_76`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-17/144*ep^-2 - 341/864*ep^-1)
        //        + cl2*sqrt3*rat(8/81*ep^-1)
        // F3/H = rat(5/48*ep^-2 + 25/96*ep^-1)
        // F5/H = rat(7/48*ep^-2 + 167/864*ep^-1)
        //        + cl2*sqrt3*rat(-8/81*ep^-1)
        // F8/H = rat(1/12*ep^-2 + 1/18*ep^-1)
        // F13/H = rat(-1/8*ep^-2 - 1/12*ep^-1)
        // F14/H = rat(-1/12*ep^-2 - 1/18*ep^-1)
        // F1..F2, F4, F6..F7, F9..F12, F15..F17 = 0
        // sum(Fi)/H = rat(1/144*ep^-2 - 7/288*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );

        let a = amp.graphs[77].renormalization_part(&settings).unwrap();
        // d78: RQFT `ghost_nnlo_77`.
        // H = p1.p1*gs^6*ca^3
        // F0/H = rat(-17/144*ep^-2 - 341/864*ep^-1)
        //        + cl2*sqrt3*rat(8/81*ep^-1)
        // F3/H = rat(5/48*ep^-2 + 25/96*ep^-1)
        // F5/H = rat(7/48*ep^-2 + 167/864*ep^-1)
        //        + cl2*sqrt3*rat(-8/81*ep^-1)
        // F8/H = rat(1/12*ep^-2 + 1/18*ep^-1)
        // F13/H = rat(-1/8*ep^-2 - 1/12*ep^-1)
        // F14/H = rat(-1/12*ep^-2 - 1/18*ep^-1)
        // F1..F2, F4, F6..F7, F9..F12, F15..F17 = 0
        // sum(Fi)/H = rat(1/144*ep^-2 - 7/288*ep^-1)
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"1/24*ca^2*dot(P(0),P(0),mink(4))*gs^4*ε^(-1)"
        );
    }
}
