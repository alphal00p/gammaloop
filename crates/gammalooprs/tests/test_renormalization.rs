use gammalooprs::{
    dot,
    graph::{Graph, parse::IntoGraph},
    initialisation::test_initialise,
    model::Model,
    processes::{Amplitude, AmplitudeGraph},
    utils::{GS, load_generic_model, symbolica_ext::LogPrint},
    uv::{
        RenormalizationPart, UVgenerationSettings,
        settings::{AlphaLoopSettings, MATADSettings, VakintSettings},
    },
};
use idenso::{
    color::{CS, ColorSimplifier},
    dirac::GammaSimplifier,
    shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip},
};
use spenso::shadowing::symbolica_utils::AtomCoreExt;
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    parse, parse_lit,
};

fn finite_part_uv_settings() -> UVgenerationSettings {
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
        only_integrated: true,
        pole_part: true,
        vakint: VakintSettings {
            normalization: undoing_normalization.to_plain_string(),
            ..Default::default()
        },
        ..Default::default()
    }
}
pub fn align_to_rqft(atom: &Atom, model: &Model) -> Atom {
    (model
        .apply_parameter_replacement_rules(
            &model.apply_coupling_replacement_rules(&-atom.simplify_color().expand()),
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
#[test]
fn finite_part_quark_lo() {
    // p1.p1*gs^2*cf*rat(ep^-1)
    // cf = 4/3
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
        .renormalization_part(&finite_part_uv_settings())
        .unwrap();

    println!("ren part: {:>}", a.log_print(Some(80)));
    insta::assert_snapshot!(
        align_to_rqft(&a,&model)
        .to_bare_ordered_string(),@"CF*dot(P(0,mink(4)),P(0,mink(4)))*gs^2*ε^(-1)"
    );
    // -1 * target
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

    let settings = finite_part_uv_settings();

    let new_settings = UVgenerationSettings {
        use_legacy: false,
        cached: false,
        ..settings.clone()
    };

    let model = load_generic_model("sm");

    #[derive(Debug)]
    struct KernelStatsSnapshot {
        cached_kernel_hits: usize,
        uncached_kernel_hits: usize,
        forest_size: usize,
    }

    impl std::fmt::Display for KernelStatsSnapshot {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(
                f,
                "cached_kernel_hits={}, uncached_kernel_hits={}, forest_size={}",
                self.cached_kernel_hits, self.uncached_kernel_hits, self.forest_size
            )
        }
    }

    fn assert_new_paths_match_legacy(
        amp: &mut AmplitudeGraph,
        a: RenormalizationPart,
        new_settings: &UVgenerationSettings,
    ) -> KernelStatsSnapshot {
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
        let uncached_kernel_hits = new_part.stats.kernel_hits;
        assert_eq!(
            normalize(&new_part.expression),
            normalize(&a.expression),
            "New renormalization for graph {} gives:\n{}\n vs old\n{}",
            amp.graph.name,
            new_part.log_print(Some(120)),
            a.log_print(Some(120))
        );

        let mut cached_settings = new_settings.clone();
        cached_settings.cached = true;

        let new_cached_part = amp.renormalization_part(&cached_settings).unwrap();
        let cached_kernel_hits = new_cached_part.stats.kernel_hits;
        assert_eq!(
            normalize(&new_cached_part.expression),
            normalize(&a.expression),
            "New cached renormalization for graph {} gives:\n{}\n vs old\n{}",
            amp.graph.name,
            new_cached_part.log_print(Some(120)),
            a.log_print(Some(120))
        );

        let cached_stats = new_cached_part.stats;
        let forest_size = cached_stats.forest_node_count;
        assert!(
            cached_kernel_hits < uncached_kernel_hits,
            "Cached path for graph {} did not reduce kernel hits: cached={}, uncached={}",
            amp.graph.name,
            cached_kernel_hits,
            uncached_kernel_hits
        );

        KernelStatsSnapshot {
            cached_kernel_hits,
            uncached_kernel_hits,
            forest_size,
        }
    }

    let a = amp.graphs[0].renormalization_part(&settings).unwrap();
    //p1.p1*i_*gs^4*CA^2*rat( - 3/16*ep^-2 + 5/32*ep^-1)
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-3𝑖/16+5𝑖/32*ε)*CA^2*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-2)");
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[0], a, &new_settings);
    insta::assert_snapshot!(
        stats.to_string(),
        @"cached_kernel_hits=5, uncached_kernel_hits=7, forest_size=6"
    );

    //p1.p1*i_*gs^4*CA^2*rat( - 1/16*ep^-2 + 1/32*ep^-1)
    let a = amp.graphs[1].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-1𝑖/32*ε+1𝑖/16)*CA^2*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-2)"
    ); //-1 * target
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[1], a, &new_settings);
    insta::assert_snapshot!(
        stats.to_string(),
        @"cached_kernel_hits=5, uncached_kernel_hits=7, forest_size=6"
    );

    //p1.p1*i_*gs^4*CA^2*rat( - 1/8*ep^-2 + 1/16*ep^-1)
    let a = amp.graphs[2].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-1𝑖/16*ε+1𝑖/8)*CA^2*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-2)"
    ); //-1 * target
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[2], a, &new_settings);
    insta::assert_snapshot!(
        stats.to_string(),
        @"cached_kernel_hits=3, uncached_kernel_hits=4, forest_size=4"
    );

    //p1.p1*i_*gs^4*CA*nf*rat(1/4*ep^-2 - 5/24*ep^-1)
    let a = amp.graphs[3].renormalization_part(&settings).unwrap();
    // Sign-flipped relative to the main snapshot after schoonschip-aware dot normalization.
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-21𝑖/8*ε+3𝑖/4)*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-2)"
    );
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[3], a, &new_settings);
    insta::assert_snapshot!(
        stats.to_string(),
        @"cached_kernel_hits=3, uncached_kernel_hits=4, forest_size=4"
    );

    //p1.p1*i_*gs^4*CA^2*rat( - 5/8*ep^-2 + 35/48*ep^-1)
    let a = amp.graphs[4].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"(37𝑖/48*ε+5𝑖/8)*CA^2*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-2)"
    ); //-1/2 * target
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[4], a, &new_settings);
    insta::assert_snapshot!(
        stats.to_string(),
        @"cached_kernel_hits=3, uncached_kernel_hits=4, forest_size=4"
    );

    //p1.p1*i_*gs^4*CA^2*rat(1/24*ep^-1)
    let a = amp.graphs[5].renormalization_part(&settings).unwrap();
    insta::assert_snapshot!(
       align_to_rqft(&a,&model).to_bare_ordered_string(),@"-5𝑖/24*CA^2*dot(P(0,mink(4)),P(0,mink(4)))*gs^4*ε^(-1)"
    );
    let stats = assert_new_paths_match_legacy(&mut amp.graphs[5], a, &new_settings);
    insta::assert_snapshot!(
        stats.to_string(),
        @"cached_kernel_hits=3, uncached_kernel_hits=4, forest_size=4"
    );
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
        .renormalization_part(&finite_part_uv_settings())
        .unwrap();

    println!("ren part: {:>}", a);
    //p1.p1*gs^2*CA*rat(1/2*ep^-1)
    // Sign-flipped relative to the main snapshot after schoonschip-aware dot normalization.
    insta::assert_snapshot!(
        align_to_rqft(&a,&model)
        .to_bare_ordered_string(),@"-1/2*CA*dot(P(0,mink(4)),P(0,mink(4)))*gs^2*ε^(-1)"
    );
}

mod failing {
    use super::*;

    fn ghost_3loop_settings() -> UVgenerationSettings {
        UVgenerationSettings {
            softct: false,
            only_integrated: true,
            pole_part: true,
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
    #[ignore = "reproduces the generated renormalization sum before alignment"]
    fn finite_part_ghost_3loop_renormalization_sum_origin_mwe() {
        test_initialise().unwrap();

        let model = load_generic_model("sm");
        let graphs: Vec<Graph> = Graph::from_path(
            concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/resources/graphs/uv_tests/rqft_ghG_3l.dot"
            ),
            &model,
        )
        .unwrap();

        let mut amp = Amplitude::from_graph_list("bub", graphs).unwrap();
        assert_eq!(amp.graphs[0].graph.name, "d1");
        let _ = amp.graphs[0]
            .renormalization_part(&ghost_3loop_settings())
            .unwrap();
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
        //p1.p1*gs^6*CA^3*rat( - 3/8*ep^-2 + 29/32*ep^-1)
        insta::assert_snapshot!(amp.graphs[0].graph.name,@"d1");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-1456*ε^2+-448/3*ε^2*𝜋^2+-896+4160/3*ε)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        //p1.p1*gs^6*CA^3*rat( - 1/16*ep^-2 + 5/192*ep^-1)
        let a = amp.graphs[1].renormalization_part(&settings).unwrap();
        insta::assert_snapshot!(amp.graphs[1].graph.name,@"d2");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"((-32+40/3*ε)*1/64*CA^3+-3/16*f(coad(8,hedge(1)),coad(8,hedge(10)),coad(8,hedge(7)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(8)))*f(coad(8,hedge(10)),coad(8,hedge(12)),coad(8,vertex(4,1)))*f(coad(8,hedge(12)),coad(8,hedge(3)),coad(8,hedge(5)))*f(coad(8,hedge(14)),coad(8,hedge(5)),coad(8,hedge(7)))*f(coad(8,hedge(14)),coad(8,hedge(8)),coad(8,vertex(4,1)))+3/16*f(coad(8,hedge(1)),coad(8,hedge(10)),coad(8,hedge(6)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(8)))*f(coad(8,hedge(10)),coad(8,hedge(12)),coad(8,vertex(4,1)))*f(coad(8,hedge(12)),coad(8,hedge(3)),coad(8,hedge(5)))*f(coad(8,hedge(14)),coad(8,hedge(5)),coad(8,hedge(6)))*f(coad(8,hedge(14)),coad(8,hedge(8)),coad(8,vertex(4,1))))*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-2)"
        );

        //p1.p1*gs^6*CA^3*rat(9/128*ep^-3 - 39/256*ep^-2 + 9/128*ep^-1)
        let a = amp.graphs[2].renormalization_part(&settings).unwrap();
        insta::assert_snapshot!(amp.graphs[2].graph.name,@"d3");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-52*ε+24+24*ε^2)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        //p1.p1*gs^6*CA^3*rat(9/128*ep^-3 - 39/256*ep^-2 + 27/128*ep^-1)
        let a = amp.graphs[3].renormalization_part(&settings).unwrap();
        insta::assert_snapshot!(amp.graphs[3].graph.name,@"d4");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-48+56*ε)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-2)"
        );

        //p1.p1*i_*gs^4*CA^2*rat( - 5/8*ep^-2 + 35/48*ep^-1)
        let a = amp.graphs[4].renormalization_part(&settings).unwrap();
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-24+-72*ε^2+52*ε)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        ); //-1/2 * target

        //p1.p1*i_*gs^4*CA^2*rat(1/24*ep^-1)
        let a = amp.graphs[5].renormalization_part(&settings).unwrap();
        insta::assert_snapshot!(amp.graphs[5].graph.name,@"d6");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-10712/3*ε+-5576*ε^2+120*ε^2*𝜋^2+3968/3*cl2*sqrt(3)*ε^2+640)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        let a = amp.graphs[6].renormalization_part(&settings).unwrap();
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-10712/3*ε+-5576*ε^2+-80+3968/3*cl2*sqrt(3)*ε^2)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        let a = amp.graphs[7].renormalization_part(&settings).unwrap();
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(-288+264*ε^2+336*ε)*1/64*CA^3*dot(P(0,mink(4)),P(0,mink(4)))*gs^6*ε^(-3)"
        );

        let a = amp.graphs[8].renormalization_part(&settings).unwrap();
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"(208*z3*ε+60*ε+64)*-1/64*dot(P(0,mink(4)),P(0,mink(4)))*f(coad(8,hedge(1)),coad(8,hedge(15)),coad(8,hedge(9)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(5)))*f(coad(8,hedge(11)),coad(8,hedge(13)),coad(8,hedge(9)))*f(coad(8,hedge(11)),coad(8,hedge(17)),coad(8,hedge(5)))*f(coad(8,hedge(13)),coad(8,hedge(3)),coad(8,hedge(7)))*f(coad(8,hedge(15)),coad(8,hedge(17)),coad(8,hedge(7)))*gs^6*ε^(-2)"
        );

        let a = amp.graphs[9].renormalization_part(&settings).unwrap();
        insta::assert_snapshot!(amp.graphs[4].graph.name,@"d5");
        insta::assert_snapshot!(
           align_to_rqft(&a,&model).to_bare_ordered_string(),@"((8*ε+8/3)*1/64*f(coad(8,hedge(1)),coad(8,hedge(11)),coad(8,hedge(15)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(5)))*f(coad(8,hedge(11)),coad(8,hedge(13)),coad(8,hedge(9)))*f(coad(8,hedge(13)),coad(8,hedge(5)),coad(8,hedge(7)))*f(coad(8,hedge(15)),coad(8,hedge(17)),coad(8,hedge(7)))+-1/16*f(coad(8,hedge(1)),coad(8,hedge(11)),coad(8,hedge(15)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(4)))*f(coad(8,hedge(11)),coad(8,hedge(13)),coad(8,hedge(9)))*f(coad(8,hedge(13)),coad(8,hedge(4)),coad(8,hedge(7)))*f(coad(8,hedge(15)),coad(8,hedge(17)),coad(8,hedge(7)))+1/16*f(coad(8,hedge(1)),coad(8,hedge(10)),coad(8,hedge(14)))*f(coad(8,hedge(1)),coad(8,hedge(3)),coad(8,hedge(5)))*f(coad(8,hedge(10)),coad(8,hedge(13)),coad(8,hedge(9)))*f(coad(8,hedge(13)),coad(8,hedge(5)),coad(8,hedge(7)))*f(coad(8,hedge(14)),coad(8,hedge(17)),coad(8,hedge(7))))*dot(P(0,mink(4)),P(0,mink(4)))*f(coad(8,hedge(17)),coad(8,hedge(3)),coad(8,hedge(9)))*gs^6*ε^(-2)"
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
}
