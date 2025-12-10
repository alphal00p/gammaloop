use bincode::de;
use color_eyre::Result;

mod test_utils;
use gammaloop_api::{commands::save, state::SyncSettings};
use gammalooprs::feyngen::diagram_generator::evaluate_sign_origin;
use gammalooprs::processes::ProcessCollection;
use gammalooprs::{feyngen::diagram_generator::evaluate_overall_factor, status_info};
use log::info;
use symbolica::{
    atom::{Atom, AtomCore},
    printer::CanonicalOrderingSettings,
};
use test_utils::{clean_test, get_test_cli, get_tests_workspace_path};
use tracing::debug;

use crate::test_utils::CLIState;
use insta::assert_snapshot;

fn count_graphs_in_processes(cli: &CLIState) -> (usize, Atom) {
    assert_eq!(cli.state.process_list.processes.len(), 1);
    let process = &cli.state.process_list.processes[0];
    let integrand_names = process.collection.get_integrand_names();
    assert_eq!(integrand_names.len(), 1);
    let integrand_name = integrand_names[0];

    let graphs = match &process.collection {
        ProcessCollection::Amplitudes(amps) => {
            let amp = amps.get(integrand_name).unwrap();
            amp.graphs.iter().map(|g| &g.graph).collect::<Vec<_>>()
        }
        ProcessCollection::CrossSections(xss) => {
            let xs = xss.get(integrand_name).unwrap();
            xs.supergraphs.iter().map(|g| &g.graph).collect::<Vec<_>>()
        }
    };
    let n_graphs = graphs.len();
    let mut overall_factor_sum = Atom::Zero;
    for g in graphs {
        overall_factor_sum += &g.overall_factor;
    }
    // overall_factor_sum = evaluate_overall_factor(overall_factor_sum.as_view());
    (n_graphs, overall_factor_sum)
}

fn split_before_flags(s: &str) -> (&str, &str) {
    match s.find("--") {
        Some(i) => {
            let (left, right) = s.split_at(i);
            (left.trim_end(), right) // X, Y
        }
        None => (s, ""), // no flags present
    }
}

fn generate_graphs_and_count(
    cli: &mut CLIState,
    generation_type: &str,
    process: &str,
    save: bool,
) -> Result<(usize, Atom)> {
    let (process_input, process_options) = split_before_flags(process);
    let cmd = format!(
        "generate {} {} {} --clear-existing-processes --only-diagrams",
        generation_type, process_input, process_options
    );
    debug!("Running diagram generation command: {}", cmd);
    cli.run_command(&cmd)?;
    if save {
        cli.run_command("save dot")?;
    }
    Ok(count_graphs_in_processes(cli))
}

fn feyngen_str(
    cli: &mut CLIState,
    generation_type: &str,
    process: &str,
    save: bool,
) -> Result<String> {
    let (n_graphs, overall_factor_sum) =
        generate_graphs_and_count(cli, generation_type, process, save)?;
    Ok(format!(
        "{} | {} = {}",
        n_graphs,
        evaluate_sign_origin(overall_factor_sum.as_view()).to_canonically_ordered_string(
            CanonicalOrderingSettings {
                include_namespace: false,
                include_attributes: false,
                ..Default::default()
            }
        ),
        evaluate_overall_factor(overall_factor_sum.as_view()).to_canonically_ordered_string(
            CanonicalOrderingSettings {
                include_namespace: false,
                include_attributes: false,
                ..Default::default()
            }
        )
    ))
}

#[test]
fn simple_epem_ddx_generation() -> Result<()> {
    let mut cli = get_test_cli(
        Some("sm_load.toml".into()),
        get_tests_workspace_path().join("epem_ddx_generation"),
        Some("epem_ddx_generation".to_string()),
        true,
    )?;

    cli.run_command("generate amp e+ e- > d d~")?;

    Ok(())
}

#[test]
fn from_symbolica() -> Result<()> {
    let mut cli = get_test_cli(
        Some("epem_ddx.toml".into()),
        get_tests_workspace_path().join("epem_ddx"),
        Some("epem_ddx".to_string()),
        true,
    )?;

    let ProcessCollection::CrossSections(x) = &cli.state.process_list.processes[0].collection
    else {
        panic!("Expected cross-section process");
    };
    let xs = x.values().next().unwrap();

    let a = xs.supergraphs.iter().next().unwrap();

    insta::assert_snapshot!(a.graph.debug_dot(),@r#"
    digraph GL0{
      num = "1";
      overall_factor = "(AutG(1))^(-1)*AntiFermionSpinSumSign(-1)*ExtFerm(-1)*IntFerm(-1)";
      overall_factor_evaluated = "-1";
      projector = "u(1,spenso::bis(4,hedge(3)))*ubar(1,spenso::bis(4,hedge(2)))*v(0,spenso::bis(4,hedge(0)))*vbar(0,spenso::bis(4,hedge(1)))";
      0[dod=0 int_id=V_71 num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(4)),spenso::dind(spenso::cof(3,hedge(6))))*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,hedge(4)),spenso::mink(4,hedge(8)))"];
      1[dod=0 int_id=V_71 num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(7)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(7)),spenso::mink(4,hedge(10)))"];
      2[dod=0 int_id=V_98 num="UFO::GC_3*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(0)),spenso::mink(4,hedge(11)))"];
      3[dod=0 int_id=V_98 num="UFO::GC_3*spenso::gamma(spenso::bis(4,hedge(1)),spenso::bis(4,hedge(3)),spenso::mink(4,hedge(9)))"];

      ext0	 [style=invis];
      ext0	-> 3:0	 [id=0 dir=back sink=0 dod=-2 is_cut=0 is_dummy=false lmb_rep="P(0,a___)" name=e0 num="1" particle="e+" pin="x:@-left,y:@edgee0"];
      ext1	 [style=invis];
      ext1	-> 3:1	 [id=1 sink=1 dod=-2 is_cut=1 is_dummy=false lmb_rep="P(1,a___)" name=e1 num="1" particle="e-" pin="x:@-left,y:@edgee1"];
      1:10	-> 2:11	 [id=2 dir=none source=2 sink=2  dod=-2 is_dummy=false lmb_rep="P(0,a___)+P(1,a___)" name=e2 num="-1*spenso::g(spenso::mink(4,hedge(10)),spenso::mink(4,hedge(11)))" particle="a"];
      0:8	-> 3:9	 [id=3 dir=none source=2 sink=2  dod=-2 is_dummy=false lmb_rep="-1*P(0,a___)+-1*P(1,a___)" name=e3 num="-1*spenso::g(spenso::mink(4,hedge(8)),spenso::mink(4,hedge(9)))" particle="a"];
      0:4	-> 1:5	 [id=4 dir=back source=1 sink=0  dod=-1 is_dummy=false lmb_id=0 lmb_rep="K(0,a___)" name=e4 num="Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(4))))*spenso::gamma(spenso::bis(4,hedge(4)),spenso::bis(4,hedge(5)),spenso::mink(4,edge(4,1)))" particle="d~"];
      0:6	-> 1:7	 [id=5 source=0 sink=1  dod=-1 is_dummy=false lmb_rep="-1*K(0,a___)+P(0,a___)+P(1,a___)" name=e5 num="Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::cof(3,hedge(6)),spenso::dind(spenso::cof(3,hedge(7))))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(6)),spenso::mink(4,edge(5,1)))" particle="d"];
      ext6	 [style=invis];
      2:3	-> ext6	 [id=6 dir=back source=1 dod=-2 is_dummy=false name=e0 num="1" particle="e+" pin="x:@+right,y:@edgee0"];
      ext7	 [style=invis];
      2:2	-> ext7	 [id=7 source=0 dod=-2 is_dummy=false name=e1 num="1" particle="e-" pin="x:@+right,y:@edgee1"];
    }
    "#);

    // clean_test(&cli.cli_settings.state_folder);

    Ok(())
}

#[test]
// fn simple_epem_ddx_generation() -> Result<()> {
//     let mut cli = get_test_cli(
//         Some("sm_load.toml".into()),
//         get_tests_workspace_path().join("epem_ddx_generation"),
//         Some("epem_ddx_generation".to_string()),
//         true,
//     )?;

//     cli.run_command("generate amp e+ e- > d d~")?;

//     clean_test(&cli.cli_settings.state_folder);

//     Ok(())
// }
#[test]
fn graph_count_from_amplitude_load() -> Result<()> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join("feyn_gen_generation_test"),
        Some("feyngen".to_string()),
        true,
    )?;
    cli.run_command("import model sm.json")?;
    cli.run_command("import graphs ./tests/resources/graphs/qqx_aaa_subtracted.dot")?;
    let (n_graphs, overall_factor_sum) = count_graphs_in_processes(&cli);
    assert_snapshot!((format!("{} | {}",n_graphs,overall_factor_sum.to_canonical_string())),@"19 | 8-11𝑖");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn example_graph_count() -> Result<()> {
    let mut cli = get_test_cli( Some("sm_load.toml".into()),get_tests_workspace_path().join("feyn_gen_generation_test"),Some("feyngen".to_string()),true)?;


    // Choose the model to consider
    // cli.run_command("import model sm.json")?;

    // A first process
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ / z QED^2==4 [{{1}} QCD] --numerator-grouping no_grouping",false)?,@"1 | -1 = -1");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ | e- a d g QED^2==4 [{{2}} QCD] --numerator-grouping only_detect_zeroes",false)?,@"3 | -3 = -3");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ | e- a d g QED^2==4 [{{2}} QCD] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"2 | -1+Group(0,1,-1)+Group(1,1,-1) = -3");//good
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e+ e- > d d~ / z QED==2 [{1}] --filter-selfenergies true",false)?,@"1 | 1 = 1");//good
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e+ e- > d d~ / z QED==2 [{1}] --filter-selfenergies false",false)?,@"3 | 3 = 3");//good

    assert_snapshot!(feyngen_str(&mut cli, "amp", "e+ e- > d d~ / z | e- a d g QED^2==4 [{{1}} QCD] --compare-canonized-numerator true --number-of-samples-for-numerator-comparisons 3",false)?,@"6 | 1+Group(1,1,1)+Group(3,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),1) = -9*G^2*Nc^(-1)*TR*ee^(-2)+2+9*G^2*Nc*TR*ee^(-2)");

    // Another process, etc...
    // ...

    Ok(())
}

#[test]
#[rustfmt::skip]
fn cp_fix_from_symbolica()->Result<()>{
    let mut cli = get_test_cli(None,get_tests_workspace_path().join("feyn_gen_generation_test"),Some("feyngen".to_string()),true)?;
    cli.cli_settings.global.logfile_directive = "[from_numerator_symbolic_expression]=debug,[compare_with_scalar_rescaling]=debug,[compare_with_sign_only]=debug".to_string();
    cli.cli_settings.sync_settings().unwrap();

    // Choose the model to consider
    cli.run_command("import model sm-default.json")?;

    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --filter-zero-flow-edges false --fully-numerical-substitution-when-comparing-numerators false --compare-canonized-numerator true",false)?,@"10 | -7+Group(10,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),-1)+Group(11,1,-1)+Group(12,1,-1)+Group(5,1,-1)+Group(6,1,-1)+Group(7,1,-1)+Group(8,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),-1)+Group(9,12*G^2*ee^(-2),-1) = -12+-12*G^2*ee^(-2)+-18*G^2*Nc*TR*ee^(-2)+18*G^2*Nc^(-1)*TR*ee^(-2)");//good
    Ok(())
}

#[test]
fn single_test() -> Result<()> {
    let mut cli = get_test_cli(
        Some("sm_load.toml".into()),
        get_tests_workspace_path().join("feyn_gen_generation_test"),
        Some("feyngen".to_string()),
        true,
    )?;

    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{3}} QCD=2] --symmetrize-left-right-states true --symmetric-left-right-polarizations true --numerator-grouping no_grouping --select-graphs GL23 GL24",true)?,@"16 | -41/2");
    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{3}} QCD=2] --symmetrize-left-right-states true --symmetric-left-right-polarizations true  --select-graphs GL23 GL24",true)?,@"0 | 0 = 0");
    println!(
        "{}",
        feyngen_str(
            &mut cli,
            "xs",
            "e+ e- > t t~ h | e+ e- g t t~ h ghG ghG~ a QCD^2==2 QED^2==6 [{{3}} QCD] --select-graphs GL06 --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --symmetrize-left-right-states true --select-graphs GL06",
            true
        )?
    );

    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --numerator-grouping no_grouping --select-graphs GL03 GL04",true)?,@"2 | (AutG(1))^(-1)*2*AntiFermionSpinSumSign(1)*ExtFerm(1)*IntFerm(-1) = -2");
    Ok(())
}
#[test]
#[rustfmt::skip]
fn test_generate_sm_a_ddx() -> Result<()> {
    let mut cli = get_test_cli(Some("sm_load.toml".into()), get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    // Choose the model to consider
    // cli.run_command("import model sm.json")?;

    // Only d, g and a as particle contents
    // Test the validity of multiple final-state specifications for cross-section generation
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > {d d~, g g} [{{1}}] | d g a --numerator-grouping only_detect_zeroes",false)?,@"1 | -1 = -1");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ g [{{2}}] | d g a --numerator-grouping only_detect_zeroes",false)?,@"3 | -3 = -3");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ g [{{2}}] | d g a --symmetrize-left-right-states true",false)?,@"2 | -1+Group(0,1,-1)+Group(1,1,-1) = -3");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ z [{{2}}] | d g a --numerator-grouping only_detect_zeroes",false)?,@"0 | 0 = 0");//good


    // Full particle contents
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{1}}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true",false)?,@"1 | -1 = -1");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states false --compare-canonized-numerator true --number-of-samples-for-numerator-comparisons 0 --fully-numerical-substitution-when-comparing-numerators false",false)?,@"12 | -9+Group(10,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),-1)+Group(5,1,-1)+Group(6,1,-1)+Group(7,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),-1)+Group(8,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),-1)+Group(9,1,-1) = -12+-27*G^2*Nc*TR*ee^(-2)+27*G^2*Nc^(-1)*TR*ee^(-2)");//not same (12 vs 10)
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true --compare-canonized-numerator true --number-of-samples-for-numerator-comparisons 0 --fully-numerical-substitution-when-comparing-numerators true",false)?,@"12 | -9+Group(10,12*G^2*ee^(-2),-1)+Group(5,1,-1)+Group(6,1,-1)+Group(7,1,-1)+Group(8,12*G^2*ee^(-2),-1)+Group(9,12*G^2*ee^(-2),-1) = -12+-36*G^2*ee^(-2)");//not same (12 vs 10)
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true --compare-canonized-numerator false --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --number-of-samples-for-numerator-comparisons 3 --fully-numerical-substitution-when-comparing-numerators false",false)?,@"10 | -7+Group(10,12*G^2*ee^(-2),-1)+Group(11,1,-1)+Group(12,1,-1)+Group(5,1,-1)+Group(6,1,-1)+Group(7,1,-1)+Group(8,12*G^2*ee^(-2),-1)+Group(9,12*G^2*ee^(-2),-1) = -12+-36*G^2*ee^(-2)");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true --compare-canonized-numerator false --number-of-samples-for-numerator-comparisons 3 --fully-numerical-substitution-when-comparing-numerators true",false)?,@"12 | -9+Group(11,1,-1)+Group(12,1,-1)+Group(5,1,-1)+Group(6,1,-1)+Group(8,1,-1)+Group(9,1,-1) = -15");//not same (12 vs 10)

    // Only 1-flavour pure QCD corrections
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{1}}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true",false)?,@"1 | -1 = -1");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{2}} QCD=1] --symmetrize-left-right-states true --symmetric-left-right-polarizations true",false)?,@"2 | -1+Group(0,1,-1)+Group(1,1,-1) = -3");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{3}} QCD=2] --symmetrize-left-right-states true --symmetric-left-right-polarizations true",false)?,@"16 | -7/2+Group(0,1,-1)+Group(1,1,-1)+Group(10,1,-2)+Group(11,1,-2)+Group(12,1,-1)+Group(13,1,-1)+Group(14,1,-1)+Group(15,1,-1)+Group(16,1,-1)+Group(17,1,-1)+Group(18,1,-1)+Group(19,1,-1)+Group(23,1,-1)+Group(24,1,-1)+Group(3,1,-1)+Group(4,1,1)+Group(5,1,1)+Group(6,1,-1)+Group(8,1,-1/2)+Group(9,1,-1/2) = -41/2");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{4}} QCD=3] --symmetrize-left-right-states true --symmetric-left-right-polarizations true",true)?,@"158 | -32/3+Group(0,1,-1)+Group(1,1,-1)+Group(10,1,1)+Group(104,1,-1)+Group(105,1,-1)+Group(106,1,1)+Group(107,1,1)+Group(108,1,-1)+Group(11,1,1)+Group(111,1,1)+Group(112,1,1)+Group(115,1,1)+Group(116,1,1)+Group(118,1,-1/4)+Group(119,1,-1/4)+Group(12,1,-1)+Group(120,1,-1)+Group(121,1,-1)+Group(122,1,-1/2)+Group(123,1,-1/2)+Group(124,1,-1/2)+Group(125,1,-1/2)+Group(126,1,-1/4)+Group(127,1,-1/4)+Group(128,1,-1/6)+Group(129,1,-1/6)+Group(13,1,-1)+Group(130,1,-1)+Group(131,1,-1)+Group(132,1,-1)+Group(133,1,-1)+Group(134,1,-2)+Group(135,1,-2)+Group(136,1,1)+Group(137,2/7,1)+Group(138,2/7,1)+Group(139,1,1)+Group(14,1,-1)+Group(140,1,1)+Group(141,2/7,1)+Group(142,1,-1)+Group(143,1,-1)+Group(144,2/7,1)+Group(145,1,1)+Group(146,1,-1)+Group(147,1,-1)+Group(148,1,-1/2)+Group(149,1,-1/2)+Group(15,1,-1)+Group(150,1,1)+Group(151,1,1)+Group(152,1,1)+Group(153,1,1)+Group(154,1,1)+Group(155,1,1)+Group(156,1,-1/2)+Group(157,1,-1/2)+Group(158,1,-1)+Group(159,1,-1)+Group(16,1,-2)+Group(160,1,-1)+Group(161,1,-1)+Group(162,1,-1/2)+Group(163,1,-1/2)+Group(164,1,-1)+Group(165,1,-1)+Group(166,1,-1/2)+Group(167,1,-1/2)+Group(168,1,-1)+Group(169,1,-1)+Group(17,1,-2)+Group(172,1,-2)+Group(173,1,-2)+Group(174,1,-1)+Group(175,1,-1)+Group(176,1,-1)+Group(177,1,-1)+Group(178,1,-1)+Group(179,1,-1)+Group(18,1,-2)+Group(180,1,-2)+Group(181,1,-2)+Group(182,1,-2)+Group(183,1,-2)+Group(184,1,-2)+Group(185,1,-2)+Group(186,1,-2)+Group(187,1,-2)+Group(188,1,-2)+Group(189,1,-2)+Group(19,1,-2)+Group(190,1,-1)+Group(191,1,-1)+Group(192,1,-2)+Group(193,1,-2)+Group(194,1,1)+Group(195,1,1)+Group(196,1,1)+Group(197,1,1)+Group(198,1,1)+Group(199,1,1)+Group(2,1,-2)+Group(20,1,-2)+Group(200,1,-1/2)+Group(201,1,-1/2)+Group(202,1,-1)+Group(203,1,-1)+Group(204,1,-1)+Group(205,1,-1)+Group(206,1,-1)+Group(207,1,-1)+Group(208,1,-1)+Group(209,1,-1)+Group(21,1,-2)+Group(215,1,-1)+Group(216,1,-1)+Group(217,1,-1)+Group(218,1,-1)+Group(22,1,-1/2)+Group(221,1,-1)+Group(222,1,-1)+Group(223,1,-1)+Group(224,1,-1)+Group(227,1,-2)+Group(228,1,-2)+Group(229,1,-2)+Group(23,1,-1/2)+Group(230,1,-2)+Group(231,1,1)+Group(232,-7/2,1)+Group(233,1,-2)+Group(234,1,-2)+Group(235,1,-2)+Group(236,1,-2)+Group(237,1,-1)+Group(238,1,-1)+Group(239,1,-1)+Group(24,1,-2)+Group(240,1,-1)+Group(241,1,-1/2)+Group(242,1,-1/2)+Group(243,1,-1)+Group(244,1,-1)+Group(245,1,-1/2)+Group(246,1,-1/2)+Group(25,1,-2)+Group(251,1,-2)+Group(252,1,-2)+Group(253,1,-2)+Group(254,1,-2)+Group(259,1,-2)+Group(26,1,-1)+Group(260,1,-2)+Group(261,1,-2)+Group(262,1,-2)+Group(263,1,-1/2)+Group(264,1,-1/2)+Group(265,1,-1)+Group(266,1,-1)+Group(267,1,-1)+Group(268,1,-1)+Group(27,1,-1)+Group(273,1,-1/2)+Group(274,1,-1/2)+Group(279,1,1)+Group(28,1,-1)+Group(280,1,1)+Group(282,1,-1)+Group(283,1,-1)+Group(284,1,-1)+Group(285,1,-1)+Group(286,1,-2)+Group(287,1,-2)+Group(288,1,-2)+Group(289,1,-2)+Group(29,1,-1)+Group(290,1,-2)+Group(291,1,-2)+Group(292,1,-2)+Group(293,1,-2)+Group(294,1,-2)+Group(295,1,-2)+Group(296,1,-2)+Group(297,1,-2)+Group(298,1,-2)+Group(299,1,-2)+Group(3,1,-2)+Group(30,1,-1)+Group(300,1,-1)+Group(301,1,-1)+Group(306,1,-2)+Group(307,1,-2)+Group(308,1,-1)+Group(309,1,-1)+Group(31,1,-1)+Group(311,1,2)+Group(312,1,2)+Group(313,1,2)+Group(314,1,2)+Group(315,1,1)+Group(316,1,1)+Group(317,1,2)+Group(318,1,2)+Group(319,1,2)+Group(32,1,-2)+Group(320,1,2)+Group(321,1,-1)+Group(322,1,-1)+Group(323,1,-1)+Group(324,1,-1)+Group(325,1,-1)+Group(326,1,-1)+Group(33,1,-2)+Group(331,1,-1/2)+Group(332,1,-1/2)+Group(333,1,-1/2)+Group(334,1,-1/2)+Group(335,1,1)+Group(336,-2/7,1)+Group(337,-2/7,1)+Group(338,1,1)+Group(339,1,-1)+Group(340,1,-1)+Group(341,1,-1)+Group(342,1,-1)+Group(343,1,-1)+Group(344,1,-1)+Group(350,1,-2)+Group(351,1,-2)+Group(352,1,-2)+Group(353,1,-2)+Group(354,1,-2)+Group(355,1,-2)+Group(360,1,-1)+Group(361,1,-1)+Group(362,1,-2)+Group(363,1,-2)+Group(38,1,-1)+Group(39,1,-1)+Group(4,1,-1)+Group(40,1,-1)+Group(41,1,-1)+Group(42,1,-1)+Group(43,1,-1)+Group(5,1,-1)+Group(57,1,-1)+Group(58,1,-1)+Group(59,1,-1)+Group(6,1,-1)+Group(60,1,-1)+Group(61,1,-1)+Group(62,1,-1)+Group(63,1,-1)+Group(64,1,-1)+Group(65,1,-1)+Group(66,1,-1)+Group(68,1,-1)+Group(69,1,-1)+Group(7,1,-1)+Group(71,1,-1)+Group(72,1,1)+Group(73,1,1)+Group(74,1,-1)+Group(75,1,-1/2)+Group(76,1,-1/2)+Group(77,1,-1)+Group(78,1,-1)+Group(79,1,2)+Group(8,1,2)+Group(80,1,2)+Group(81,1,1)+Group(82,1,1)+Group(83,1,1)+Group(84,1,1)+Group(85,1,1)+Group(86,1,1)+Group(87,1,2)+Group(88,1,2)+Group(89,1,1/2)+Group(9,1,2)+Group(90,1,1/2)+Group(91,1,1)+Group(92,1,1)+Group(93,1,-1)+Group(94,1,1)+Group(95,1,1)+Group(96,1,-1)+Group(97,1,1)+Group(98,1,1)+Group(99,1,-1) = -2939/14");//less 158 vs 166 (lorentz cancellations)


    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_slow_generate_sm_a_ddx() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Full particle contents
    // Adding --symmetrize_left_right_states below would error in Python; keep as-is here.
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{3}}] --numerator-grouping only_detect_zeroes",false)?,@"1321 | -1103/2 = -1103/2");//good

    // Only 1-flavour pure QCD corrections
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{5}}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"6303 | -51683/24 = -51683/24");//good

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_sm_full_a_ddx() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm-full.json"))?;

    // Full particle contents
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{1}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"1 | -1 = -1");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"47 | -47 = -47");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"37 | -35+Group(29,1,-1)+Group(30,1,-1)+Group(31,1,-1)+Group(32,1,-1) = -39");//less 37 vs 45 due to lorentz cancellations
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true --numerator-grouping group_identical_graphs_up_to_sign",true)?,@"36 | -33+Group(29,1,-1)+Group(30,1,-1)+Group(32,1,-1)+Group(33,1,-1)+Group(35,1,-1)+Group(36,1,-1) = -39");//as above

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_slow_generate_sm_full_a_ddx() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm-full.json"))?;

    // Full particle contents
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ QED^2==4 [{{3}}] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"291 | -274+Group(324,1,-1)+Group(325,1,-1)+Group(328,1,-1)+Group(329,1,-1)+Group(330,1,-1)+Group(331,1,-1)+Group(336,1,-1)+Group(337,1,-1)+Group(338,1,-1)+Group(339,1,-1)+Group(340,1,-1)+Group(341,1,-1)+Group(346,1,-1)+Group(351,1,-1)+Group(352,1,-1)+Group(353,1,-1)+Group(354,1,-1)+Group(355,1,-1)+Group(356,1,-1)+Group(357,1,-1)+Group(358,1,-1)+Group(359,1,-1)+Group(360,1,-1)+Group(365,1,-1)+Group(366,1,-1)+Group(367,1,-1)+Group(368,1,-1)+Group(369,1,-1)+Group(370,1,-1)+Group(371,1,-1)+Group(38,1,-1)+Group(39,1,-1)+Group(40,1,-1)+Group(41,1,-1)+Group(44,1,-1)+Group(45,1,-1)+Group(46,1,-1)+Group(47,1,-1)+Group(500,1,-1)+Group(501,1,-1)+Group(502,1,-1)+Group(503,1,-1)+Group(504,1,-1)+Group(505,1,-1)+Group(506,1,-1)+Group(507,1,-1) = -320");//a lot less TODO: check that main has lorentz zeros present and here removed 291 vs 339
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{3}}] --numerator-grouping no_grouping",false)?,@"8549 | -9111/2 = -9111/2");//good

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_generate_sm_h_n_j() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Full particle contents
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g | h g b t ghg [{1}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"2 | -4 = -4");//good
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g | h g b t ghg [{2}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"30 | -34+Group(1,1,-2)+Group(13,1,-2)+Group(14,1,-2)+Group(15,1,-2)+Group(16,1,-2)+Group(2,1,-2)+Group(21,1,-2)+Group(22,1,-2)+Group(25,1,-2)+Group(26,1,-2)+Group(27,1,-2)+Group(28,1,-2)+Group(3,1,-2)+Group(33,1,-2)+Group(34,1,-2)+Group(35,1,-2)+Group(36,1,-2)+Group(4,1,-2)+Group(5,1,-2)+Group(6,1,-2) = -74");//good
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g | h g b t ghg [{1}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"2 | -4 = -4");//good
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g | h g b t ghg [{2}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"40 | -74 = -74");//good
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g g | h g b t ghg [{1}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"6 | -24 = -24");//good
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g g | h g b t ghg [{2}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"188 | -618 = -618");//good
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g g | h g b t ghg [{1}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"4 | -12+Group(1,1,-3)+Group(2,1,-3)+Group(4,1,-3)+Group(5,1,-3) = -24");//good
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g g | h g b t ghg [{2}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"144 | -480+Group(10,1,6)+Group(100,1,-3/2)+Group(101,1,-3/2)+Group(102,1,6)+Group(103,1,6)+Group(104,1,6)+Group(105,1,6)+Group(106,1,-6)+Group(107,1,-6)+Group(108,1,-3)+Group(109,1,-3)+Group(11,1,6)+Group(110,1,-3)+Group(111,1,-3)+Group(112,1,-3/2)+Group(113,1,-3/2)+Group(114,1,3)+Group(115,1,3)+Group(120,1,3)+Group(121,1,3)+Group(122,1,-3)+Group(123,1,-3)+Group(127,1,-3)+Group(128,1,-3)+Group(129,1,-3)+Group(130,1,-3)+Group(136,1,6)+Group(137,1,6)+Group(162,1,-3)+Group(163,1,-3)+Group(170,1,-3)+Group(171,1,-3)+Group(172,1,-3)+Group(173,1,-3)+Group(174,1,-3)+Group(175,1,-3)+Group(176,1,-3)+Group(177,1,-3)+Group(178,1,-3)+Group(179,1,-3)+Group(180,1,-3)+Group(181,1,-3)+Group(189,1,-3/2)+Group(190,1,-3/2)+Group(191,1,-3)+Group(192,1,-3)+Group(33,1,-3)+Group(34,1,-3)+Group(36,1,-3)+Group(37,1,-3)+Group(38,1,-3)+Group(39,1,-3)+Group(4,1,3)+Group(47,1,-3)+Group(48,1,-3)+Group(5,1,3)+Group(52,1,-3)+Group(53,1,-3)+Group(54,1,-6)+Group(55,1,-6)+Group(56,1,-3)+Group(57,1,-3)+Group(58,1,-3)+Group(59,1,-3)+Group(6,1,3)+Group(60,1,-3/2)+Group(61,1,-3/2)+Group(63,1,-3/2)+Group(64,1,-3/2)+Group(65,1,-3)+Group(66,1,-3)+Group(7,1,3)+Group(71,1,-3)+Group(72,1,-3)+Group(73,1,-3)+Group(74,1,-3)+Group(84,1,-3)+Group(85,1,-3)+Group(88,1,-3)+Group(89,1,-3)+Group(90,1,-3)+Group(91,1,-3)+Group(94,1,-3)+Group(95,1,-3)+Group(96,1,-3)+Group(97,1,-3)+Group(98,1,-3/2)+Group(99,1,-3/2) = -618");// less grouping, but same sum 132 vs 144 ?

    // Cross-section at 3-loops
    assert_snapshot!(feyngen_str(&mut cli, "xs", "h > g g | h g b t ghg [{{3}}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"8 | 8 = 8");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "h > g g | h g b t ghg [{{3}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"3 | Group(0,1,1)+Group(1,1,1)+Group(2,1,1)+Group(3,1,1)+Group(4,1,1)+Group(5,1,1)+Group(6,1,1)+Group(7,1,1) = 8");//good

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_slow_generate_sm_h_n_j_cross_section() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "xs", "h > g g g | h g b t ghg [{{4}}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"72 | 88 = 88");//good
    assert_snapshot!(feyngen_str(&mut cli, "xs", "h > g g g | h g b t ghg [{{4}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"26 | 4+Group(10,1,1)+Group(11,1,1)+Group(12,1,1/2)+Group(13,1,1/2)+Group(14,1,1/2)+Group(15,1,1/2)+Group(16,1,2)+Group(17,1,2)+Group(18,1,2)+Group(19,1,2)+Group(2,1,1)+Group(20,1,1)+Group(21,1,1)+Group(22,1,1)+Group(23,1,1)+Group(24,1,1)+Group(25,1,1)+Group(26,1,1)+Group(27,1,1)+Group(28,1,2)+Group(29,1,2)+Group(3,1,1)+Group(30,1,2)+Group(31,1,2)+Group(32,1,1)+Group(33,1,1)+Group(34,1,1)+Group(35,1,1)+Group(36,1,1)+Group(37,1,1)+Group(38,1,1)+Group(39,1,1)+Group(4,1,1)+Group(40,1,2)+Group(41,1,2)+Group(42,1,2)+Group(43,1,2)+Group(44,1,1)+Group(45,1,1)+Group(46,1,1)+Group(47,1,1)+Group(48,1,2)+Group(49,1,2)+Group(5,1,1)+Group(50,1,2)+Group(51,1,2)+Group(54,1,1)+Group(55,1,1)+Group(56,1,1)+Group(57,1,1)+Group(58,1,1)+Group(59,1,1)+Group(6,1,1)+Group(60,1,1/2)+Group(61,1,1/2)+Group(62,1,1/2)+Group(63,1,1/2)+Group(64,1,1)+Group(65,1,1)+Group(66,1,1)+Group(67,1,1)+Group(68,1,2)+Group(69,1,2)+Group(7,1,1)+Group(70,1,2)+Group(71,1,2)+Group(8,1,1)+Group(9,1,1) = 88");//good

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_slow_generate_sm_h_n_j() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g | h g b t ghg [{3}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true --numerator-grouping only_detect_zeroes",false)?,@"1194 | -1260 = -1260");//good
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g g | h g b t ghg [{3}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true --numerator-grouping only_detect_zeroes",false)?,@"7018 | -12747 = -12747");//sum is larger and there are more graphs.. this is questionable

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_epem_a_qqh() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Cross sections
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"3 | -3 = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"3 | -3 = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"30 | -30 = -30");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] --numerator-grouping only_detect_zeroes",false)?,@"594 | -261 = -261");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"21 | -21 = -21");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=2] --numerator-grouping only_detect_zeroes",false)?,@"171 | -171 = -171");

    // Enable grouping
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -1+Group(0,1,-1)+Group(1,1,-1) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -1+Group(0,1,-1)+Group(1,1,-1) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"16 | -2+Group(10,1,-1)+Group(11,1,-1)+Group(12,1,-1)+Group(13,1,-1)+Group(14,1,-1)+Group(15,1,-1)+Group(16,1,-1)+Group(17,1,-1)+Group(18,1,-1)+Group(19,1,-1)+Group(22,1,-1)+Group(23,1,-1)+Group(28,1,-1)+Group(29,1,-1)+Group(30,1,-1)+Group(31,1,-1)+Group(32,1,-1)+Group(33,1,-1)+Group(34,1,-1)+Group(35,1,-1)+Group(4,1,-1)+Group(44,1,-1)+Group(45,1,-1)+Group(48,1,-1)+Group(49,1,-1)+Group(5,1,-1)+Group(50,1,-1)+Group(51,1,-1) = -30");
    // This grouping involves 21 new cancellations between graphs, so that modifies the total overall factor
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"266 | -7+Group(100,1,-1)+Group(1000,1,1)+Group(101,1,-1)+Group(1019,1,-1)+Group(102,1,1)+Group(1020,1,-1)+Group(1021,1,-1)+Group(1022,1,-1)+Group(1023,1,-1)+Group(1024,1,-1)+Group(1025,1,-1)+Group(1026,1,-1)+Group(1027,1,-1)+Group(1028,1,-1)+Group(1029,1,-1)+Group(103,1,1)+Group(1030,1,-1)+Group(1031,1,-1)+Group(1032,1,-1)+Group(1033,1,-1)+Group(1034,1,-1)+Group(1035,1,-1)+Group(1036,1,-1)+Group(1037,1,-1)+Group(1038,1,-1)+Group(104,1,1)+Group(1041,1,-1)+Group(1042,1,-1)+Group(1043,1,-1)+Group(1044,1,-1)+Group(1045,1,-1)+Group(1046,1,-1)+Group(1049,1,-1)+Group(105,1,1)+Group(1050,1,-1)+Group(1051,1,-1)+Group(1052,1,-1)+Group(1053,1,-1)+Group(1054,1,-1)+Group(1055,1,-1)+Group(1056,1,-1)+Group(1057,1,-1)+Group(1058,1,-1)+Group(1059,1,-1)+Group(1060,1,-1)+Group(1065,1,-1)+Group(1066,1,-1)+Group(1067,1,-1)+Group(1068,1,-1)+Group(1069,1,-1)+Group(1070,1,-1)+Group(1071,1,-1)+Group(1076,1,-1)+Group(1077,1,-1)+Group(1078,1,-1)+Group(1079,1,-1)+Group(1080,1,-1)+Group(1081,1,-1)+Group(1082,1,-1)+Group(1083,1,-1)+Group(1084,1,-1)+Group(1089,1,-1)+Group(1090,1,-1)+Group(1091,1,-1)+Group(1092,1,-1)+Group(1099,1,-1)+Group(1100,1,-1)+Group(1101,1,-1)+Group(1102,1,-1)+Group(130,1,-1)+Group(131,1,-1)+Group(132,1,-1)+Group(133,1,-1)+Group(134,1,-1)+Group(135,1,-1)+Group(136,1,-1)+Group(137,1,-1)+Group(138,1,-1)+Group(139,1,-1)+Group(140,1,-1)+Group(141,1,-1)+Group(142,1,-1)+Group(143,1,-1)+Group(144,1,-1)+Group(145,1,-1)+Group(148,1,-1)+Group(149,1,-1)+Group(150,1,-1)+Group(151,1,-1)+Group(152,1,-1)+Group(153,1,-1)+Group(154,1,-1)+Group(155,1,-1)+Group(160,1,-1)+Group(161,1,-1)+Group(162,1,-1)+Group(163,1,-1)+Group(164,1,-1)+Group(165,1,-1)+Group(166,1,-1)+Group(167,1,-1)+Group(172,1,1)+Group(173,1,1)+Group(174,1,1)+Group(175,1,1)+Group(176,1,1)+Group(177,1,1)+Group(178,1,1)+Group(179,1,1)+Group(180,1,-1)+Group(181,1,-1)+Group(182,1,-1)+Group(183,1,-1)+Group(184,1,-1)+Group(185,1,-1)+Group(186,1,-1)+Group(187,1,-1)+Group(188,1,-1)+Group(189,1,-1)+Group(190,1,-1)+Group(191,1,-1)+Group(192,1,-1)+Group(193,1,-1)+Group(194,1,-1)+Group(195,1,-1)+Group(20,1,-1)+Group(204,1,1)+Group(205,1,1)+Group(206,1,1)+Group(207,1,1)+Group(208,1,1)+Group(209,1,1)+Group(21,1,-1)+Group(22,1,1)+Group(224,1,-1/2)+Group(225,1,-1/2)+Group(226,1,-1/2)+Group(227,1,-1/2)+Group(228,1,-1/2)+Group(229,1,-1/2)+Group(23,1,1)+Group(230,1,-1/2)+Group(231,1,-1/2)+Group(232,1,-1)+Group(233,1,-1)+Group(234,1,-1)+Group(235,1,-1)+Group(236,1,-1)+Group(237,1,-1)+Group(238,1,-1)+Group(239,1,-1)+Group(24,1,-1/2)+Group(240,1,-1)+Group(241,1,-1)+Group(242,1,-1)+Group(243,1,-1)+Group(244,1,-1)+Group(245,1,-1)+Group(246,1,-1)+Group(247,1,-1)+Group(248,1,-1)+Group(249,1,-1)+Group(25,1,-1/2)+Group(250,1,-1)+Group(251,1,-1)+Group(252,1,-1)+Group(253,1,-1)+Group(254,1,-1)+Group(255,1,-1)+Group(256,1,-1)+Group(257,1,-1)+Group(258,1,-1)+Group(259,1,-1)+Group(26,1,-1)+Group(27,1,-1)+Group(28,1,-1)+Group(280,1,-1)+Group(281,1,1)+Group(282,1,1)+Group(283,1,-1)+Group(284,1,-1)+Group(285,1,-1)+Group(286,1,1)+Group(287,1,1)+Group(288,1,-1/2)+Group(289,1,-1/2)+Group(29,1,-1)+Group(30,1,-1)+Group(302,1,-1)+Group(303,1,-1)+Group(304,1,-1)+Group(305,1,-1)+Group(306,1,-1)+Group(307,1,-1)+Group(308,1,-1)+Group(309,1,-1)+Group(31,1,-1)+Group(314,1,-1)+Group(315,1,-1)+Group(316,1,-1)+Group(317,1,-1)+Group(318,1,-1)+Group(319,1,-1)+Group(32,1,1)+Group(33,1,1)+Group(330,1,-1)+Group(331,1,-1)+Group(332,1,1)+Group(333,1,1)+Group(34,1,-1)+Group(342,1,-1)+Group(343,1,-1)+Group(344,1,1)+Group(345,1,1)+Group(346,1,-1)+Group(347,1,-1)+Group(348,1,1)+Group(349,1,1)+Group(35,1,-1)+Group(36,1,-1)+Group(37,1,-1)+Group(38,1,-1)+Group(386,1,-1)+Group(387,1,-1)+Group(388,1,-1)+Group(389,1,-1)+Group(39,1,-1)+Group(390,1,1)+Group(391,1,1)+Group(392,1,1)+Group(393,1,1)+Group(402,1,-1)+Group(403,1,-1)+Group(406,1,-1)+Group(407,1,-1)+Group(408,1,-1)+Group(409,1,-1)+Group(410,1,-1)+Group(411,1,-1)+Group(412,1,1)+Group(413,1,1)+Group(414,1,1)+Group(415,1,1)+Group(423,1,-1)+Group(424,1,-1)+Group(425,1,-1)+Group(426,1,-1)+Group(427,1,-1)+Group(428,1,-1)+Group(447,1,-1/2)+Group(448,1,-1/2)+Group(449,1,-1/2)+Group(450,1,-1/2)+Group(46,1,-1)+Group(47,1,-1)+Group(471,1,-1)+Group(472,1,-1)+Group(473,1,-1)+Group(474,1,-1)+Group(479,1,-1)+Group(48,1,-1)+Group(480,1,-1)+Group(481,1,-1)+Group(482,1,-1)+Group(49,1,-1)+Group(50,1,-1)+Group(51,1,-1)+Group(523,1,-1)+Group(524,1,-1)+Group(525,1,-1)+Group(526,1,-1)+Group(527,1,-1)+Group(528,1,-1)+Group(529,1,-1)+Group(530,1,-1)+Group(539,1,-1)+Group(540,1,-1)+Group(541,1,-1)+Group(542,1,-1)+Group(543,1,-1)+Group(544,1,-1)+Group(545,1,-1)+Group(546,1,-1)+Group(563,1,-1)+Group(564,1,-1)+Group(565,1,-1)+Group(566,1,-1)+Group(567,1,-1)+Group(568,1,-1)+Group(569,1,-1)+Group(570,1,-1)+Group(571,1,-1)+Group(572,1,-1)+Group(573,1,-1)+Group(574,1,-1)+Group(575,1,-1)+Group(576,1,-1)+Group(577,1,-1)+Group(578,1,-1)+Group(579,1,-1)+Group(580,1,-1)+Group(581,1,-1)+Group(582,1,-1)+Group(583,1,-1)+Group(584,1,-1)+Group(585,1,-1)+Group(586,1,-1)+Group(587,1,1)+Group(588,1,1)+Group(589,1,1)+Group(590,1,1)+Group(591,1,1)+Group(592,1,1)+Group(593,1,1)+Group(594,1,1)+Group(595,1,-1)+Group(596,1,-1)+Group(597,1,-1)+Group(598,1,-1)+Group(599,1,1)+Group(600,1,1)+Group(601,1,-1)+Group(602,1,-1)+Group(603,1,1)+Group(604,1,1)+Group(605,1,1)+Group(606,1,1)+Group(607,1,-1)+Group(608,1,-1)+Group(609,1,1)+Group(610,1,1)+Group(611,1,-1)+Group(612,1,-1)+Group(613,1,-1)+Group(614,1,-1)+Group(641,1,-1)+Group(642,1,-1)+Group(643,1,-1)+Group(644,1,-1)+Group(649,1,-1)+Group(650,1,-1)+Group(651,1,-1)+Group(652,1,-1)+Group(657,1,-1)+Group(658,1,-1)+Group(659,1,-1)+Group(660,1,-1)+Group(665,1,-1)+Group(666,1,-1)+Group(667,1,-1)+Group(668,1,-1)+Group(669,1,-1)+Group(670,1,-1)+Group(671,1,-1)+Group(672,1,-1)+Group(673,1,-1)+Group(674,1,-1)+Group(675,1,-1)+Group(676,1,-1)+Group(685,1,-1)+Group(686,1,-1)+Group(687,1,-1)+Group(688,1,-1)+Group(693,1,-1)+Group(694,1,-1)+Group(695,1,-1)+Group(696,1,-1)+Group(697,1,-1)+Group(698,1,-1)+Group(699,1,-1)+Group(700,1,-1)+Group(709,1,-1)+Group(710,1,-1)+Group(711,1,-1)+Group(712,1,-1)+Group(713,1,-1)+Group(714,1,-1)+Group(715,1,-1)+Group(716,1,-1)+Group(717,1,-1)+Group(718,1,-1)+Group(719,1,-1)+Group(720,1,-1)+Group(733,1,-1)+Group(734,1,-1)+Group(735,1,-1)+Group(736,1,-1)+Group(737,1,-1)+Group(738,1,-1)+Group(739,1,-1)+Group(740,1,-1)+Group(749,1,-1)+Group(750,1,-1)+Group(751,1,-1)+Group(752,1,-1)+Group(757,1,-1)+Group(758,1,-1)+Group(759,1,-1)+Group(760,1,-1)+Group(761,1,-1)+Group(762,1,-1)+Group(763,1,-1)+Group(764,1,-1)+Group(769,1,-1)+Group(770,1,-1)+Group(771,1,-1)+Group(772,1,-1)+Group(775,1,1)+Group(776,1,1)+Group(777,1,1)+Group(778,1,1)+Group(78,1,-1)+Group(79,1,-1)+Group(791,1,1)+Group(792,1,1)+Group(80,1,-1)+Group(809,1,-1/2)+Group(81,1,-1)+Group(810,1,-1/2)+Group(815,1,-1/2)+Group(816,1,-1/2)+Group(817,1,-1/2)+Group(818,1,-1/2)+Group(819,1,-1)+Group(82,1,1)+Group(820,1,-1)+Group(821,1,1)+Group(822,1,1)+Group(823,1,-1/2)+Group(824,1,-1/2)+Group(83,1,1)+Group(835,1,-1)+Group(836,1,-1)+Group(837,1,-1)+Group(838,1,-1)+Group(839,1,1)+Group(84,1,1)+Group(840,1,1)+Group(841,1,1)+Group(842,1,1)+Group(85,1,1)+Group(859,1,-1)+Group(860,1,-1)+Group(873,1,-1)+Group(874,1,-1)+Group(875,1,-1)+Group(876,1,-1)+Group(881,1,-1)+Group(882,1,-1)+Group(883,1,-1)+Group(884,1,-1)+Group(889,1,-1)+Group(890,1,-1)+Group(891,1,-1)+Group(892,1,-1)+Group(909,1,-1/2)+Group(910,1,-1/2)+Group(911,1,-1/2)+Group(912,1,-1/2)+Group(919,1,-1)+Group(920,1,-1)+Group(921,1,-1)+Group(922,1,-1)+Group(923,1,-1)+Group(924,1,-1)+Group(925,1,-1)+Group(926,1,-1)+Group(935,1,-1)+Group(936,1,-1)+Group(941,1,1)+Group(942,1,-1)+Group(943,1,-1)+Group(944,1,1)+Group(949,1,1)+Group(950,1,1)+Group(951,1,1)+Group(952,1,1)+Group(953,1,-1)+Group(954,1,-1)+Group(955,1,-1)+Group(956,1,-1)+Group(957,1,-1)+Group(958,1,-1)+Group(959,1,-1)+Group(960,1,-1)+Group(969,1,-1)+Group(970,1,-1)+Group(971,1,1)+Group(972,1,1)+Group(975,1,1)+Group(976,1,1)+Group(977,1,1)+Group(978,1,1)+Group(979,1,-1)+Group(98,1,-1)+Group(980,1,-1)+Group(981,1,-1)+Group(982,1,-1)+Group(983,1,-1)+Group(984,1,-1)+Group(986,1,-1)+Group(989,1,-1)+Group(99,1,-1)+Group(991,1,-1)+Group(992,1,-1)+Group(995,1,1)+Group(996,1,1)+Group(997,1,-1)+Group(998,1,-1)+Group(999,1,1) = -327");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"20 | -19+Group(0,1,-1)+Group(1,1,-1) = -21");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=1] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"157 | -143+Group(100,1,-1)+Group(101,1,-1)+Group(110,1,-1)+Group(111,1,-1)+Group(114,1,-1)+Group(115,1,-1)+Group(116,1,-1)+Group(117,1,-1)+Group(12,1,-1)+Group(13,1,-1)+Group(18,1,-1)+Group(19,1,-1)+Group(20,1,-1)+Group(21,1,-1)+Group(22,1,-1)+Group(23,1,-1)+Group(24,1,-1)+Group(25,1,-1)+Group(26,1,-1)+Group(27,1,-1)+Group(88,1,-1)+Group(89,1,-1)+Group(94,1,-1)+Group(95,1,-1)+Group(96,1,-1)+Group(97,1,-1)+Group(98,1,-1)+Group(99,1,-1) = -171");

    // Enable left-right-symmetry
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -1+Group(0,1,-1)+Group(1,1,-1) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -1+Group(0,1,-1)+Group(1,1,-1) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"11 | Group(0,1,-1)+Group(1,1,-1)+Group(10,1,-2)+Group(11,1,-2)+Group(12,1,-1)+Group(13,1,-1)+Group(18,1,-1)+Group(19,1,-1)+Group(2,1,-2)+Group(20,1,-1)+Group(21,1,-1)+Group(22,1,-1)+Group(23,1,-1)+Group(26,1,-1)+Group(27,1,-1)+Group(3,1,-2)+Group(30,1,-2)+Group(31,1,-2)+Group(6,1,-2)+Group(7,1,-2)+Group(8,1,-1)+Group(9,1,-1) = -30");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"166 | -2+Group(0,1,-1)+Group(1,1,-1)+Group(10,1,1)+Group(100,1,-2)+Group(101,1,-2)+Group(102,1,-1)+Group(103,1,-1)+Group(104,1,-1)+Group(105,1,-1)+Group(106,1,-1)+Group(107,1,-1)+Group(108,1,-2)+Group(109,1,-2)+Group(11,1,-1)+Group(110,1,-2)+Group(111,1,-2)+Group(112,1,-2)+Group(113,1,-2)+Group(114,1,1)+Group(115,1,1)+Group(116,1,-1)+Group(117,1,-1)+Group(12,1,-1)+Group(128,1,-2)+Group(129,1,-2)+Group(13,1,-1)+Group(130,1,-2)+Group(131,1,-2)+Group(132,1,-2)+Group(133,1,-2)+Group(134,1,-1)+Group(135,1,-1)+Group(136,1,-1)+Group(137,1,-1)+Group(138,1,-1)+Group(139,1,-1)+Group(14,1,-1)+Group(140,1,1)+Group(141,1,1)+Group(142,1,1)+Group(143,1,1)+Group(144,1,-2)+Group(145,1,-2)+Group(146,1,-1)+Group(147,1,-1)+Group(148,1,2)+Group(149,1,2)+Group(15,1,-1)+Group(150,1,-2)+Group(151,1,-2)+Group(152,1,2)+Group(153,1,-2)+Group(154,1,-2)+Group(155,1,2)+Group(156,1,-2)+Group(157,1,-2)+Group(158,1,-2)+Group(159,1,-2)+Group(16,1,-2)+Group(160,1,-1)+Group(161,1,-1)+Group(162,1,-1/2)+Group(163,1,-1/2)+Group(164,1,1)+Group(165,1,1)+Group(166,1,-1)+Group(167,1,-1)+Group(168,1,1)+Group(169,1,-1)+Group(17,1,-2)+Group(170,1,-1)+Group(171,1,1)+Group(172,1,-1)+Group(173,1,-1)+Group(174,1,-1)+Group(175,1,-1)+Group(18,1,-1)+Group(188,1,-1)+Group(189,1,-1)+Group(19,1,-1)+Group(190,1,-1/2)+Group(191,1,-1/2)+Group(192,1,1)+Group(193,1,1)+Group(194,1,-1)+Group(195,1,-1)+Group(196,1,-1/2)+Group(197,1,-1/2)+Group(198,1,1)+Group(199,1,1)+Group(2,1,-1/2)+Group(20,1,2)+Group(200,1,-1)+Group(201,1,-1)+Group(206,1,-2)+Group(207,1,-2)+Group(208,1,-1)+Group(209,1,-1)+Group(21,1,2)+Group(210,1,-1)+Group(211,1,-1)+Group(212,1,-1)+Group(213,1,-1)+Group(214,1,-1)+Group(215,1,-1)+Group(217,1,-1)+Group(22,1,-2)+Group(220,1,-1)+Group(222,1,-1)+Group(223,1,-1)+Group(224,1,-1)+Group(229,1,-1)+Group(23,1,-2)+Group(230,1,-1)+Group(231,1,-1)+Group(232,1,-2)+Group(233,1,-2)+Group(238,1,-2)+Group(239,1,-2)+Group(24,1,2)+Group(248,1,-2)+Group(249,1,-2)+Group(25,1,-2)+Group(250,1,-2)+Group(251,1,-2)+Group(252,1,-2)+Group(253,1,-2)+Group(254,1,-1)+Group(255,1,-1)+Group(26,1,-2)+Group(260,1,1)+Group(261,1,1)+Group(262,1,1)+Group(263,1,1)+Group(264,1,-1)+Group(265,1,-1)+Group(266,1,-1/2)+Group(267,1,-1/2)+Group(268,1,1)+Group(269,1,1)+Group(27,1,2)+Group(278,1,-2)+Group(279,1,-2)+Group(28,1,-2)+Group(284,1,-1)+Group(285,1,-1)+Group(288,1,-1)+Group(289,1,-1)+Group(29,1,-2)+Group(290,1,-1)+Group(291,1,-1)+Group(298,1,-2)+Group(299,1,-2)+Group(3,1,-1/2)+Group(30,1,-2)+Group(300,1,-2)+Group(301,1,-2)+Group(302,1,-2)+Group(303,1,-2)+Group(308,1,1)+Group(309,1,1)+Group(31,1,-2)+Group(312,1,-1)+Group(313,1,-1)+Group(314,1,-1/2)+Group(315,1,-1/2)+Group(316,1,1)+Group(317,1,1)+Group(326,1,-2)+Group(327,1,-2)+Group(330,1,1)+Group(331,1,-1)+Group(332,1,1)+Group(333,1,1)+Group(334,1,-1)+Group(335,1,1)+Group(336,1,-1)+Group(337,1,-1)+Group(338,1,-1)+Group(339,1,-1)+Group(344,1,-2)+Group(345,1,-2)+Group(346,1,-2)+Group(347,1,-2)+Group(348,1,-2)+Group(349,1,-2)+Group(356,1,1)+Group(357,1,1)+Group(360,1,-2)+Group(361,1,-2)+Group(362,1,-1)+Group(363,1,-1)+Group(364,1,2)+Group(365,1,2)+Group(376,1,-2)+Group(377,1,-2)+Group(378,1,-2)+Group(379,1,-2)+Group(38,1,-2)+Group(383,1,-2)+Group(386,1,-2)+Group(388,1,-2)+Group(389,1,-2)+Group(39,1,-2)+Group(390,1,-2)+Group(391,1,-2)+Group(4,1,1)+Group(40,1,-1)+Group(400,1,-2)+Group(401,1,-2)+Group(402,1,-2)+Group(403,1,-2)+Group(404,1,-2)+Group(405,1,-2)+Group(406,1,-2)+Group(407,1,-2)+Group(41,1,-1)+Group(410,1,-2)+Group(411,1,-2)+Group(412,1,-2)+Group(413,1,-2)+Group(418,1,2)+Group(419,1,2)+Group(42,1,2)+Group(422,1,-2)+Group(423,1,-2)+Group(43,1,2)+Group(44,1,-2)+Group(442,1,-2)+Group(443,1,-2)+Group(45,1,-2)+Group(450,1,-2)+Group(451,1,-2)+Group(452,1,-2)+Group(453,1,-2)+Group(454,1,-2)+Group(455,1,-2)+Group(456,1,-2)+Group(457,1,-2)+Group(458,1,-2)+Group(459,1,-2)+Group(460,1,-2)+Group(461,1,-2)+Group(474,1,-2)+Group(475,1,-2)+Group(476,1,-2)+Group(477,1,-2)+Group(478,1,-2)+Group(479,1,-2)+Group(480,1,-2)+Group(481,1,-2)+Group(484,1,-2)+Group(485,1,-2)+Group(486,1,-2)+Group(487,1,-2)+Group(492,1,2)+Group(493,1,2)+Group(494,1,2)+Group(495,1,2)+Group(5,1,1)+Group(500,1,-2)+Group(501,1,-2)+Group(502,1,-2)+Group(503,1,-2)+Group(504,1,-2)+Group(505,1,-2)+Group(52,1,-2)+Group(520,1,-2)+Group(521,1,-2)+Group(524,1,-2)+Group(525,1,-2)+Group(526,1,-2)+Group(527,1,-2)+Group(53,1,-2)+Group(532,1,-1)+Group(533,1,-1)+Group(536,1,-1)+Group(537,1,-1)+Group(538,1,-1)+Group(539,1,-1)+Group(54,1,-2)+Group(540,1,-1)+Group(541,1,-1)+Group(546,1,-2)+Group(547,1,-2)+Group(55,1,-2)+Group(550,1,-2)+Group(551,1,-2)+Group(552,1,-2)+Group(553,1,-2)+Group(56,1,-2)+Group(565,1,-2)+Group(566,1,-2)+Group(567,1,-1)+Group(568,1,-1)+Group(57,1,-2)+Group(58,1,-2)+Group(59,1,2)+Group(6,1,-1)+Group(60,1,2)+Group(61,1,2)+Group(62,1,2)+Group(63,1,-2)+Group(64,1,-2)+Group(65,1,-2)+Group(66,1,-2)+Group(67,1,-2)+Group(68,1,-2)+Group(69,1,-2)+Group(693,1,1)+Group(694,1,1)+Group(695,1,1)+Group(696,1,1)+Group(697,1,1)+Group(698,1,1)+Group(7,1,-1)+Group(701,1,1)+Group(702,1,1)+Group(8,1,-1)+Group(80,1,-2)+Group(81,1,-2)+Group(86,1,-2)+Group(87,1,-2)+Group(88,1,-2)+Group(89,1,-2)+Group(9,1,1)+Group(92,1,2)+Group(93,1,2)+Group(94,1,-1)+Group(95,1,-1)+Group(96,1,-1/2)+Group(97,1,-1/2)+Group(98,1,1)+Group(99,1,1) = -327");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"16 | -19+Group(0,1,-1)+Group(1,1,-1) = -21");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=1] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"104 | -141+Group(0,1,-1)+Group(1,1,-1)+Group(10,1,-2)+Group(11,1,-2)+Group(12,1,-1)+Group(13,1,-1)+Group(18,1,-1)+Group(19,1,-1)+Group(2,1,-2)+Group(20,1,-1)+Group(21,1,-1)+Group(22,1,-1)+Group(23,1,-1)+Group(26,1,-1)+Group(27,1,-1)+Group(3,1,-2)+Group(30,1,-2)+Group(31,1,-2)+Group(6,1,-2)+Group(7,1,-2)+Group(8,1,-1)+Group(9,1,-1) = -171");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_a_qqh() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Cross sections
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d a ghg g QED^2==2 [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"3 | -3 = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"3 | -3 = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{3}} QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"30 | -30 = -30");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{4}} QCD=2] --numerator-grouping only_detect_zeroes",false)?,@"594 | -261 = -261");

    // Enable left-right-symmetry and grouping
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d a ghg g QED^2==2 [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -1+Group(0,1,-1)+Group(1,1,-1) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -1+Group(0,1,-1)+Group(1,1,-1) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{3}} QCD=1] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"11 | Group(0,1,-1)+Group(1,1,-1)+Group(10,1,-1)+Group(11,1,-1)+Group(12,1,-2)+Group(13,1,-2)+Group(14,1,-2)+Group(15,1,-2)+Group(16,1,-2)+Group(17,1,-2)+Group(18,1,-1)+Group(19,1,-1)+Group(2,1,-2)+Group(20,1,-1)+Group(21,1,-1)+Group(26,1,-1)+Group(27,1,-1)+Group(28,1,-1)+Group(29,1,-1)+Group(3,1,-2)+Group(4,1,-1)+Group(5,1,-1) = -30");
    // Additional cancellations modify the total overall factor
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{4}} QCD=2] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"166 | -2+Group(0,1,-1)+Group(1,1,1)+Group(10,1,-2)+Group(102,1,-1)+Group(103,1,-1)+Group(104,1,-1)+Group(105,1,-1)+Group(106,1,-1)+Group(107,1,-1)+Group(108,1,-1)+Group(109,1,-1)+Group(11,1,-2)+Group(110,1,-2)+Group(111,1,-2)+Group(112,1,-2)+Group(113,1,-2)+Group(114,1,-1)+Group(115,1,-1)+Group(116,1,2)+Group(117,1,2)+Group(118,1,-1)+Group(119,1,-1)+Group(120,1,1)+Group(121,1,1)+Group(122,1,-1/2)+Group(123,1,-1/2)+Group(124,1,-1)+Group(125,1,-1)+Group(126,1,-2)+Group(127,1,-2)+Group(132,1,-2)+Group(133,1,-2)+Group(134,1,-1)+Group(135,1,-1)+Group(16,1,-2)+Group(166,1,1)+Group(167,1,1)+Group(168,1,-1)+Group(169,1,-1)+Group(17,1,-2)+Group(170,1,1)+Group(171,1,1)+Group(172,1,-2)+Group(173,1,-2)+Group(174,1,2)+Group(175,1,2)+Group(176,1,2)+Group(177,1,2)+Group(18,1,-1)+Group(184,1,-1)+Group(185,1,-1)+Group(186,1,-1)+Group(187,1,-1)+Group(188,1,-2)+Group(189,1,-2)+Group(19,1,1)+Group(190,1,-2)+Group(191,1,-2)+Group(192,1,-2)+Group(193,1,-2)+Group(194,1,-1)+Group(195,1,-1)+Group(196,1,-1/2)+Group(197,1,-1/2)+Group(198,1,1)+Group(199,1,1)+Group(2,1,1)+Group(20,1,1)+Group(200,1,-2)+Group(201,1,-2)+Group(202,1,-2)+Group(203,1,-2)+Group(204,1,-1)+Group(205,1,-1)+Group(206,1,2)+Group(207,1,2)+Group(21,1,-1)+Group(212,1,-2)+Group(213,1,-2)+Group(214,1,-2)+Group(215,1,-2)+Group(220,1,-2)+Group(221,1,-2)+Group(222,1,-2)+Group(223,1,-2)+Group(224,1,-2)+Group(225,1,-2)+Group(226,1,-1)+Group(227,1,-1)+Group(228,1,2)+Group(229,1,2)+Group(234,1,-2)+Group(235,1,-2)+Group(240,1,-2)+Group(241,1,-2)+Group(242,1,-2)+Group(243,1,-2)+Group(244,1,-2)+Group(245,1,-2)+Group(246,1,-2)+Group(247,1,-2)+Group(252,1,-2)+Group(253,1,-2)+Group(254,1,-2)+Group(255,1,-2)+Group(256,1,-2)+Group(257,1,-2)+Group(26,1,1)+Group(262,1,-2)+Group(263,1,-2)+Group(264,1,-2)+Group(265,1,-2)+Group(27,1,1)+Group(28,1,1)+Group(282,1,-2)+Group(283,1,-2)+Group(284,1,-2)+Group(285,1,-2)+Group(286,1,-2)+Group(287,1,-2)+Group(29,1,1)+Group(3,1,-1)+Group(30,1,-1)+Group(302,1,-2)+Group(303,1,-2)+Group(304,1,-2)+Group(305,1,-2)+Group(306,1,-2)+Group(307,1,-2)+Group(308,1,-2)+Group(309,1,-2)+Group(31,1,-1)+Group(310,1,-2)+Group(311,1,-2)+Group(312,1,-2)+Group(313,1,-2)+Group(314,1,-2)+Group(315,1,-2)+Group(316,1,-2)+Group(317,1,-2)+Group(318,1,-2)+Group(319,1,-2)+Group(32,1,-2)+Group(320,1,-2)+Group(321,1,-2)+Group(322,1,-2)+Group(323,1,-2)+Group(324,1,-2)+Group(325,1,-2)+Group(326,1,-2)+Group(327,1,2)+Group(328,1,2)+Group(329,1,2)+Group(33,1,-2)+Group(330,1,2)+Group(331,1,-2)+Group(332,1,-2)+Group(333,1,-2)+Group(34,1,-1)+Group(340,1,1)+Group(341,1,1)+Group(342,1,2)+Group(343,1,-2)+Group(344,1,-2)+Group(345,1,2)+Group(346,1,-2)+Group(347,1,-2)+Group(348,1,2)+Group(349,1,2)+Group(35,1,-1)+Group(350,1,2)+Group(351,1,2)+Group(356,1,-2)+Group(357,1,-2)+Group(358,1,-2)+Group(359,1,-2)+Group(36,1,-1)+Group(360,1,-2)+Group(361,1,-2)+Group(362,1,-2)+Group(363,1,-2)+Group(364,1,-2)+Group(365,1,-2)+Group(366,1,-2)+Group(367,1,-2)+Group(368,1,-2)+Group(369,1,-2)+Group(37,1,-1)+Group(370,1,-2)+Group(371,1,-2)+Group(372,1,-2)+Group(373,1,-2)+Group(374,1,-2)+Group(375,1,-2)+Group(38,1,-2)+Group(380,1,-2)+Group(381,1,-2)+Group(382,1,2)+Group(383,1,2)+Group(384,1,-1)+Group(385,1,-1)+Group(386,1,1)+Group(387,1,1)+Group(388,1,1)+Group(389,1,1)+Group(39,1,-2)+Group(390,1,1)+Group(391,1,1)+Group(392,1,-1)+Group(393,1,-1)+Group(394,1,-1)+Group(395,1,-1)+Group(396,1,-1/2)+Group(397,1,-1/2)+Group(398,1,1)+Group(399,1,1)+Group(4,1,2)+Group(40,1,-2)+Group(404,1,-1)+Group(405,1,-1)+Group(406,1,-1)+Group(407,1,-1)+Group(408,1,-1)+Group(409,1,-1)+Group(41,1,-2)+Group(410,1,-2)+Group(411,1,-2)+Group(412,1,-2)+Group(413,1,-2)+Group(414,1,-2)+Group(415,1,-2)+Group(416,1,-2)+Group(417,1,-2)+Group(418,1,-2)+Group(419,1,-2)+Group(42,1,-2)+Group(420,1,-2)+Group(421,1,-2)+Group(422,1,-2)+Group(423,1,-2)+Group(424,1,1)+Group(425,1,1)+Group(426,1,-1)+Group(427,1,-1)+Group(428,1,-1)+Group(429,1,-1)+Group(43,1,-2)+Group(430,1,-1/2)+Group(431,1,-1/2)+Group(432,1,1)+Group(433,1,1)+Group(434,1,-2)+Group(435,1,-2)+Group(436,1,-1)+Group(437,1,-1)+Group(438,1,-1)+Group(439,1,-1)+Group(44,1,-2)+Group(444,1,1)+Group(445,1,1)+Group(447,1,-2)+Group(448,1,-2)+Group(449,1,-1)+Group(45,1,-2)+Group(450,1,-1)+Group(451,1,-1)+Group(452,1,-1)+Group(46,1,-1)+Group(468,1,-2)+Group(469,1,-2)+Group(47,1,-1)+Group(470,1,-2)+Group(471,1,-2)+Group(472,1,-2)+Group(473,1,-2)+Group(474,1,-1)+Group(475,1,-1)+Group(476,1,-1)+Group(477,1,-1)+Group(478,1,-1)+Group(479,1,-1)+Group(48,1,-1/2)+Group(485,1,-1)+Group(486,1,-1)+Group(487,1,-1)+Group(488,1,-1)+Group(49,1,-1/2)+Group(493,1,-1/2)+Group(494,1,-1/2)+Group(495,1,-1/2)+Group(496,1,-1/2)+Group(5,1,2)+Group(50,1,1)+Group(501,1,1)+Group(502,1,1)+Group(503,1,1)+Group(504,1,1)+Group(505,1,-2)+Group(506,1,-2)+Group(51,1,1)+Group(513,1,1)+Group(514,1,1)+Group(515,1,1)+Group(516,1,1)+Group(56,1,-1)+Group(57,1,-1)+Group(58,1,-1)+Group(59,1,-1)+Group(6,1,1)+Group(60,1,-1)+Group(65,1,-1)+Group(68,1,-1)+Group(69,1,-1)+Group(7,1,1)+Group(72,1,-1)+Group(73,1,-1)+Group(74,1,-1)+Group(75,1,-1)+Group(8,1,-1)+Group(80,1,-2)+Group(81,1,-2)+Group(82,1,-2)+Group(83,1,-2)+Group(84,1,-2)+Group(85,1,-2)+Group(9,1,-1)+Group(90,1,-2)+Group(91,1,-2)+Group(92,1,-1)+Group(93,1,-1)+Group(94,1,-1)+Group(95,1,-1)+Group(96,1,-1)+Group(97,1,-1) = -327");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | b h a ghg g QED^2==4 [{{4}} QCD=2] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"151 | -3+Group(0,1,-1)+Group(1,1,1)+Group(10,1,-2)+Group(100,1,-1)+Group(101,1,-1)+Group(102,1,-1)+Group(103,1,-1)+Group(104,1,-2)+Group(105,1,-2)+Group(106,1,-2)+Group(107,1,-2)+Group(108,1,-1)+Group(109,1,-1)+Group(11,1,-2)+Group(110,1,-1)+Group(111,1,-1)+Group(112,1,-1/2)+Group(113,1,-1/2)+Group(114,1,-1)+Group(115,1,-1)+Group(116,1,-2)+Group(117,1,-2)+Group(122,1,-2)+Group(123,1,-2)+Group(124,1,-1)+Group(125,1,-1)+Group(148,1,1)+Group(149,1,1)+Group(150,1,-1)+Group(151,1,-1)+Group(152,1,1)+Group(153,1,1)+Group(154,1,-2)+Group(155,1,-2)+Group(156,1,2)+Group(157,1,2)+Group(158,1,2)+Group(159,1,2)+Group(16,1,-2)+Group(166,1,-1)+Group(167,1,-1)+Group(168,1,-1)+Group(169,1,-1)+Group(17,1,-2)+Group(170,1,-2)+Group(171,1,-2)+Group(172,1,-2)+Group(173,1,-2)+Group(174,1,-2)+Group(175,1,-2)+Group(176,1,-1)+Group(177,1,-1)+Group(178,1,-1/2)+Group(179,1,-1/2)+Group(18,1,-1)+Group(180,1,-2)+Group(181,1,-2)+Group(182,1,-2)+Group(183,1,-2)+Group(184,1,-1)+Group(185,1,-1)+Group(186,1,-2)+Group(187,1,-2)+Group(188,1,-2)+Group(189,1,-2)+Group(19,1,1)+Group(194,1,-2)+Group(195,1,-2)+Group(196,1,-2)+Group(197,1,-2)+Group(198,1,-2)+Group(199,1,-2)+Group(2,1,1)+Group(20,1,1)+Group(200,1,-1)+Group(201,1,-1)+Group(202,1,-2)+Group(203,1,-2)+Group(208,1,-2)+Group(209,1,-2)+Group(21,1,-1)+Group(210,1,-2)+Group(211,1,-2)+Group(212,1,-2)+Group(213,1,-2)+Group(214,1,-2)+Group(215,1,-2)+Group(220,1,-2)+Group(221,1,-2)+Group(222,1,-2)+Group(223,1,-2)+Group(224,1,-2)+Group(225,1,-2)+Group(230,1,-2)+Group(231,1,-2)+Group(232,1,-2)+Group(233,1,-2)+Group(250,1,-2)+Group(251,1,-2)+Group(252,1,-2)+Group(253,1,-2)+Group(254,1,-2)+Group(255,1,-2)+Group(26,1,1)+Group(27,1,1)+Group(270,1,-2)+Group(271,1,-2)+Group(272,1,-2)+Group(273,1,-2)+Group(274,1,-2)+Group(275,1,-2)+Group(276,1,-2)+Group(277,1,-2)+Group(278,1,-2)+Group(279,1,-2)+Group(28,1,1)+Group(280,1,-2)+Group(281,1,-2)+Group(282,1,-2)+Group(283,1,-2)+Group(284,1,-2)+Group(285,1,-2)+Group(286,1,-2)+Group(287,1,-2)+Group(288,1,-2)+Group(289,1,-2)+Group(29,1,1)+Group(290,1,-2)+Group(291,1,-2)+Group(292,1,-2)+Group(293,1,-2)+Group(294,1,-2)+Group(295,1,2)+Group(296,1,2)+Group(297,1,2)+Group(298,1,2)+Group(299,1,-2)+Group(3,1,-1)+Group(30,1,-1)+Group(300,1,-2)+Group(301,1,-2)+Group(308,1,1)+Group(309,1,1)+Group(31,1,-1)+Group(310,1,2)+Group(311,1,-2)+Group(312,1,-2)+Group(313,1,2)+Group(314,1,-2)+Group(315,1,-2)+Group(316,1,2)+Group(317,1,2)+Group(318,1,2)+Group(319,1,2)+Group(32,1,-2)+Group(324,1,-2)+Group(325,1,-2)+Group(326,1,-2)+Group(327,1,-2)+Group(328,1,-2)+Group(329,1,-2)+Group(33,1,-2)+Group(330,1,-2)+Group(331,1,-2)+Group(332,1,-2)+Group(333,1,-2)+Group(334,1,-2)+Group(335,1,-2)+Group(336,1,-2)+Group(337,1,-2)+Group(338,1,-2)+Group(339,1,-2)+Group(34,1,-1)+Group(340,1,-2)+Group(341,1,-2)+Group(342,1,-2)+Group(343,1,-2)+Group(348,1,-2)+Group(349,1,-2)+Group(35,1,-1)+Group(350,1,-1)+Group(351,1,-1)+Group(352,1,-1)+Group(353,1,-1)+Group(354,1,-1)+Group(355,1,-1)+Group(356,1,-1/2)+Group(357,1,-1/2)+Group(358,1,-1)+Group(359,1,-1)+Group(36,1,-1)+Group(360,1,-1)+Group(361,1,-1)+Group(362,1,-1)+Group(363,1,-1)+Group(364,1,-2)+Group(365,1,-2)+Group(366,1,-2)+Group(367,1,-2)+Group(368,1,-2)+Group(369,1,-2)+Group(37,1,-1)+Group(370,1,-2)+Group(371,1,-2)+Group(372,1,-2)+Group(373,1,-2)+Group(374,1,-2)+Group(375,1,-2)+Group(376,1,-2)+Group(377,1,-2)+Group(378,1,1)+Group(379,1,1)+Group(38,1,-2)+Group(380,1,-1)+Group(381,1,-1)+Group(382,1,-1)+Group(383,1,-1)+Group(384,1,-1/2)+Group(385,1,-1/2)+Group(386,1,-2)+Group(387,1,-2)+Group(388,1,-1)+Group(389,1,-1)+Group(39,1,-2)+Group(390,1,-1)+Group(391,1,-1)+Group(396,1,-2)+Group(397,1,-2)+Group(398,1,-1)+Group(399,1,-1)+Group(4,1,2)+Group(40,1,-2)+Group(400,1,-1)+Group(401,1,-1)+Group(41,1,-2)+Group(417,1,-2)+Group(418,1,-2)+Group(419,1,-2)+Group(42,1,-2)+Group(420,1,-2)+Group(421,1,-2)+Group(422,1,-2)+Group(423,1,-1)+Group(424,1,-1)+Group(425,1,-1)+Group(426,1,-1)+Group(427,1,-1)+Group(428,1,-1)+Group(43,1,-2)+Group(434,1,-1)+Group(435,1,-1)+Group(436,1,-1)+Group(437,1,-1)+Group(44,1,-2)+Group(442,1,-1/2)+Group(443,1,-1/2)+Group(444,1,-1/2)+Group(445,1,-1/2)+Group(446,1,-2)+Group(447,1,-2)+Group(45,1,-2)+Group(454,1,1)+Group(455,1,1)+Group(456,1,1)+Group(457,1,1)+Group(46,1,-1)+Group(47,1,-1)+Group(48,1,-1/2)+Group(49,1,-1/2)+Group(5,1,2)+Group(50,1,-1)+Group(51,1,-1)+Group(52,1,-1)+Group(53,1,-1)+Group(54,1,-1)+Group(59,1,-1)+Group(6,1,1)+Group(62,1,-1)+Group(63,1,-1)+Group(66,1,-1)+Group(67,1,-1)+Group(68,1,-1)+Group(69,1,-1)+Group(7,1,1)+Group(74,1,-2)+Group(75,1,-2)+Group(76,1,-2)+Group(77,1,-2)+Group(78,1,-2)+Group(79,1,-2)+Group(8,1,-1)+Group(84,1,-2)+Group(85,1,-2)+Group(86,1,-1)+Group(87,1,-1)+Group(88,1,-1)+Group(89,1,-1)+Group(9,1,-1)+Group(90,1,-1)+Group(91,1,-1)+Group(96,1,-1)+Group(97,1,-1)+Group(98,1,-1)+Group(99,1,-1) = -366");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_aa_ttx_amplitude() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{0} QCD=0] --numerator-grouping only_detect_zeroes",false)?,@"2 | 2 = 2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{1} QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"8 | 8 = 8");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{2} QCD=2] --numerator-grouping only_detect_zeroes",false)?,@"152 | 28 = 28");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{3} QCD=3] --numerator-grouping only_detect_zeroes",false)?,@"3608 | -104/3 = -104/3");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{0} QCD=0] --symmetrize-initial-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"1 | 2 = 2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{1} QCD=1] --symmetrize-initial-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"4 | 8 = 8");
    // Six additional cancellations imply that the total overall factor is different
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{2} QCD=2] --symmetrize-initial-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"58 | 68+Group(56,1,-2)+Group(57,1,-2)+Group(90,1,-2)+Group(91,1,-2) = 60");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{3} QCD=3] --symmetrize-initial-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"1201 | 2548/3+Group(101,1,-2)+Group(102,1,-2)+Group(1024,1,-2)+Group(1025,-7/2,-2)+Group(1028,1,-2)+Group(1029,-2/7,-2)+Group(1054,1,-2)+Group(1055,7/2,-2)+Group(1071,1,-2)+Group(1072,1,-2)+Group(1080,1,-2)+Group(1081,1,-2)+Group(1091,1,-2)+Group(1092,-2/7,-2)+Group(1093,1,-2)+Group(1094,2/7,-2)+Group(1124,1,-2)+Group(1125,-7/2,-2)+Group(1126,1,-2)+Group(1127,7/2,-2)+Group(1164,1,-2)+Group(1165,-2/7,-2)+Group(1166,1,-2)+Group(1167,1,-2)+Group(1179,1,-2)+Group(1180,-7/2,-2)+Group(1202,1,-2)+Group(1203,1,-2)+Group(1204,1,-2)+Group(1205,1,-2)+Group(1209,1,-2)+Group(1210,-7/2,-2)+Group(1211,1,-2)+Group(1212,1,-2)+Group(1224,1,-2)+Group(1225,-2/7,-2)+Group(1247,1,-2)+Group(1248,1,-2)+Group(1249,1,-2)+Group(1250,1,-2)+Group(1253,1,-2)+Group(1254,1,-2)+Group(128,1,2)+Group(129,1,2)+Group(1329,1,-2)+Group(1330,1,-2)+Group(1331,1,-2)+Group(1332,1,-2)+Group(1333,1,-1)+Group(1334,1,-1)+Group(1384,1,-2)+Group(1385,1,-2)+Group(1393,1,-2)+Group(1394,1,-2)+Group(1397,1,-2)+Group(1398,1,-2)+Group(1399,1,-2)+Group(1400,1,-2)+Group(1437,1,-2)+Group(1438,1,-2)+Group(1439,1,-2)+Group(1440,1,-2)+Group(1441,1,-1)+Group(1442,1,-1)+Group(1539,1,2)+Group(1540,1,2)+Group(1541,1,2)+Group(1542,1,2)+Group(1543,1,2)+Group(1544,1,2)+Group(155,1,2)+Group(156,1,2)+Group(161,1,-2)+Group(163,1,-2)+Group(1686,1,2)+Group(1687,1,2)+Group(1707,1,2)+Group(1708,1,2)+Group(1744,1,-2)+Group(1745,1,-2)+Group(1761,1,2)+Group(1762,1,2)+Group(1767,1,2)+Group(1768,1,2)+Group(1785,1,-2)+Group(1786,1,-2)+Group(1797,1,-2)+Group(1798,7/2,-2)+Group(1805,1,-2)+Group(1806,2/7,-2)+Group(1813,1,-2)+Group(1814,1,-2)+Group(1821,1,2)+Group(1822,1,2)+Group(1827,1,2)+Group(1828,1,2)+Group(1851,1,-2)+Group(1852,1,-2)+Group(1855,1,-2)+Group(1856,1,-2)+Group(1871,1,-2)+Group(1872,7/2,-2)+Group(1887,1,-2)+Group(1888,1,-2)+Group(1898,1,-2)+Group(1899,1,-2)+Group(1902,1,-2)+Group(1903,1,-2)+Group(1912,1,-2)+Group(1913,1,-2)+Group(1914,1,-1)+Group(1915,1,-1)+Group(1916,1,-2)+Group(1917,1,-2)+Group(1940,1,-2)+Group(1941,1,-2)+Group(1942,1,-1)+Group(1943,1,-1)+Group(1944,1,-2)+Group(1945,1,-2)+Group(1960,1,-2)+Group(1961,1,-2)+Group(2000,1,-2)+Group(2001,1,-2)+Group(2002,1,-2)+Group(2003,1,-2)+Group(2013,1,-2)+Group(2014,2/7,-2)+Group(2015,1,-2)+Group(2016,7/2,-2)+Group(2017,1,-2)+Group(2018,2/7,-2)+Group(2021,1,-2)+Group(2022,7/2,-2)+Group(2023,1,-2)+Group(2024,7/2,-2)+Group(2025,1,-2)+Group(2026,2/7,-2)+Group(2027,1,-2)+Group(2028,2/7,-2)+Group(2033,1,-2)+Group(2034,7/2,-2)+Group(2062,1,-2)+Group(2063,-2/7,-2)+Group(2064,1,-2)+Group(2065,-7/2,-2)+Group(2066,1,-2)+Group(2067,7/2,-2)+Group(2075,1,-2)+Group(2076,-7/2,-2)+Group(2077,1,-2)+Group(2078,-2/7,-2)+Group(2080,1,-2)+Group(2081,-7/2,-2)+Group(2082,1,-2)+Group(2083,-2/7,-2)+Group(2086,1,-2)+Group(2087,-2/7,-2)+Group(2088,1,-2)+Group(2089,2/7,-2)+Group(2090,1,-2)+Group(2091,-7/2,-2)+Group(2092,1,-2)+Group(2093,7/2,-2)+Group(2094,1,-2)+Group(2095,-2/7,-2)+Group(2096,1,-2)+Group(2097,1,-2)+Group(2098,1,-2)+Group(2099,-7/2,-2)+Group(2100,1,-2)+Group(2101,-7/2,-2)+Group(2102,1,-2)+Group(2103,1,-2)+Group(2104,1,-2)+Group(2105,-2/7,-2)+Group(2106,1,-2)+Group(2107,1,-2)+Group(2121,1,-2)+Group(2122,1,-2)+Group(2130,1,-2)+Group(2131,1,-2)+Group(2155,1,-2)+Group(2156,1,-2)+Group(2161,1,-2)+Group(2162,1,-2)+Group(217,1,-2)+Group(218,2/7,-2)+Group(2183,1,-2)+Group(2184,1,-2)+Group(2189,1,-2)+Group(2190,1,-2)+Group(2199,1,-2)+Group(220,1,-2)+Group(2200,1,-2)+Group(2201,1,-2)+Group(2202,1,-2)+Group(2203,1,-2)+Group(2204,1,-2)+Group(2205,1,-2)+Group(2206,1,-2)+Group(221,7/2,-2)+Group(2215,1,-2)+Group(2216,1,-2)+Group(222,1,-2)+Group(2228,1,-2)+Group(2229,1,-2)+Group(224,2/7,-2)+Group(2270,1,-2)+Group(2271,1,-2)+Group(2272,1,-2)+Group(2273,1,-2)+Group(232,1,-2)+Group(233,7/2,-2)+Group(235,1,-2)+Group(236,7/2,-2)+Group(239,1,-2)+Group(241,2/7,-2)+Group(242,1,-2)+Group(244,2/7,-2)+Group(251,1,-2)+Group(252,7/2,-2)+Group(448,1,-2)+Group(449,1,-2)+Group(450,1,-2)+Group(451,1,-2)+Group(452,1,-2)+Group(453,1,-2)+Group(454,1,-2)+Group(455,1,-2)+Group(470,1,-2)+Group(471,1,-2)+Group(534,1,-2)+Group(535,-2/7,-2)+Group(556,1,-2)+Group(557,1,-2)+Group(571,1,-2)+Group(572,-7/2,-2)+Group(575,1,-2)+Group(576,7/2,-2)+Group(613,1,-2)+Group(614,1,-2)+Group(625,1,-2)+Group(626,1,-2)+Group(637,1,-2)+Group(638,1,-2)+Group(643,1,-2)+Group(644,1,-2)+Group(665,1,-2)+Group(666,1,-2)+Group(671,1,-2)+Group(672,1,-2)+Group(735,1,2)+Group(736,1,2)+Group(910,1,-2)+Group(911,-7/2,-2)+Group(914,1,-2)+Group(915,-2/7,-2)+Group(941,1,-2)+Group(942,7/2,-2)+Group(958,1,-2)+Group(959,2/7,-2)+Group(972,1,-2)+Group(973,1,-2) = 1504/3");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_aa_ttx_cross_section() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --numerator-grouping only_detect_zeroes",false)?,@"4 | -4 = -4");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"40 | -40 = -40");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] --numerator-grouping only_detect_zeroes",false)?,@"874 | -266 = -266");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --symmetrize-initial-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | Group(0,1,-1)+Group(1,1,-1)+Group(2,1,-1)+Group(3,1,-1) = -4");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --symmetrize-initial-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"14 | Group(0,1,-1)+Group(1,1,-1)+Group(10,1,-1)+Group(11,1,-1)+Group(12,1,-2)+Group(13,1,-2)+Group(14,1,-2)+Group(15,1,-2)+Group(16,1,-1)+Group(17,1,-1)+Group(18,1,-1)+Group(19,1,-1)+Group(2,1,-1)+Group(24,1,-1)+Group(25,1,-1)+Group(26,1,-1)+Group(27,1,-1)+Group(28,1,-2)+Group(29,1,-2)+Group(3,1,-1)+Group(30,1,-2)+Group(31,1,-2)+Group(4,1,-2)+Group(5,1,-2)+Group(6,1,-2)+Group(7,1,-2)+Group(8,1,-1)+Group(9,1,-1) = -40");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] --symmetrize-initial-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"223 | Group(0,1,-1)+Group(1,1,1)+Group(10,1,2)+Group(100,1,-2)+Group(101,1,-2)+Group(102,1,-2)+Group(103,1,-2)+Group(104,1,-2)+Group(105,1,-2)+Group(11,1,-2)+Group(110,1,-2)+Group(111,1,-2)+Group(112,1,-2)+Group(113,1,-2)+Group(114,1,-2)+Group(115,1,-2)+Group(116,1,-2)+Group(117,1,-2)+Group(118,1,-2)+Group(119,1,-2)+Group(12,1,-2)+Group(120,1,-2)+Group(121,1,-2)+Group(122,1,-2)+Group(123,1,-2)+Group(124,1,-2)+Group(125,1,-2)+Group(126,1,-2)+Group(127,1,-2)+Group(128,1,-2)+Group(129,1,-2)+Group(13,1,2)+Group(130,1,-2)+Group(131,1,-2)+Group(132,1,-2)+Group(133,1,-2)+Group(134,1,-2)+Group(135,1,-2)+Group(14,1,2)+Group(140,1,-2)+Group(141,1,-2)+Group(142,1,-2)+Group(143,1,-2)+Group(144,1,-2)+Group(145,1,-2)+Group(146,1,-2)+Group(147,1,-2)+Group(148,1,-2)+Group(149,1,-2)+Group(15,1,-2)+Group(150,1,-2)+Group(151,1,-2)+Group(152,1,-2)+Group(153,1,-2)+Group(154,1,-2)+Group(155,1,-2)+Group(16,1,-1)+Group(164,1,-1)+Group(165,1,-1)+Group(166,1,-1)+Group(167,1,-1)+Group(168,1,-2)+Group(169,1,-2)+Group(17,1,1)+Group(170,1,-2)+Group(171,1,-2)+Group(172,1,-1)+Group(173,1,-1)+Group(174,1,-1)+Group(175,1,-1)+Group(176,1,-2)+Group(177,1,-2)+Group(178,1,-2)+Group(179,1,-2)+Group(18,1,1)+Group(180,1,-2)+Group(181,1,-2)+Group(186,1,-1)+Group(187,1,-1)+Group(188,1,-1)+Group(189,1,-1)+Group(19,1,-1)+Group(190,1,-2)+Group(191,1,-2)+Group(192,1,-2)+Group(193,1,-2)+Group(194,1,-2)+Group(195,1,-2)+Group(196,1,-2)+Group(197,1,-2)+Group(2,1,1)+Group(20,1,-1)+Group(202,1,-1)+Group(203,1,-1)+Group(204,1,-1)+Group(205,1,-1)+Group(206,1,-2)+Group(207,1,-2)+Group(208,1,-2)+Group(209,1,-2)+Group(21,1,1)+Group(210,1,-2)+Group(211,1,-2)+Group(22,1,1)+Group(220,1,-1)+Group(221,1,-1)+Group(222,1,-1)+Group(223,1,-1)+Group(224,1,-2)+Group(225,1,-2)+Group(226,1,-2)+Group(227,1,-2)+Group(228,1,-1)+Group(229,1,-1)+Group(23,1,-1)+Group(230,1,-1)+Group(231,1,-1)+Group(232,1,-1)+Group(233,1,-1)+Group(234,1,-1)+Group(235,1,-1)+Group(236,1,-2)+Group(237,1,-2)+Group(238,1,-2)+Group(239,1,-2)+Group(24,1,2)+Group(240,1,-1)+Group(241,1,-1)+Group(242,1,-1)+Group(243,1,-1)+Group(246,1,-2)+Group(247,1,-2)+Group(25,1,2)+Group(254,1,-2)+Group(255,1,-2)+Group(256,1,-1)+Group(257,1,-1)+Group(258,1,-1)+Group(259,1,-1)+Group(26,1,2)+Group(260,1,-2)+Group(261,1,-2)+Group(262,1,-1)+Group(263,1,-1)+Group(264,1,-1)+Group(265,1,-1)+Group(266,1,1)+Group(267,1,1)+Group(268,1,1)+Group(269,1,1)+Group(27,1,2)+Group(270,1,1)+Group(271,1,1)+Group(272,1,1)+Group(273,1,1)+Group(274,1,-1)+Group(275,1,-1)+Group(276,1,-1)+Group(277,1,-1)+Group(278,1,-1)+Group(279,1,-1)+Group(28,1,1)+Group(280,1,-1)+Group(281,1,-1)+Group(286,1,-2)+Group(287,1,-2)+Group(288,1,-2)+Group(289,1,-2)+Group(29,1,1)+Group(294,1,-2)+Group(295,1,-2)+Group(296,1,-2)+Group(297,1,-2)+Group(3,1,-1)+Group(30,1,1)+Group(306,1,-1)+Group(307,1,-1)+Group(308,1,-1)+Group(309,1,-1)+Group(31,1,1)+Group(310,1,-2)+Group(311,1,-2)+Group(312,1,-2)+Group(313,1,-2)+Group(314,1,-2)+Group(315,1,-2)+Group(316,1,-1)+Group(317,1,-1)+Group(318,1,-1)+Group(319,1,-1)+Group(324,1,-2)+Group(325,1,-2)+Group(326,1,-2)+Group(327,1,-2)+Group(328,1,-2)+Group(329,1,-2)+Group(330,1,-2)+Group(331,1,-2)+Group(332,1,-2)+Group(333,1,-2)+Group(334,1,1)+Group(335,1,1)+Group(336,1,1)+Group(337,1,1)+Group(338,1,-2)+Group(339,1,-2)+Group(340,1,-2)+Group(341,1,-2)+Group(342,1,-2)+Group(343,1,-2)+Group(344,1,-2)+Group(345,1,-2)+Group(346,1,1)+Group(347,1,1)+Group(348,1,1)+Group(349,1,1)+Group(350,1,-2)+Group(351,1,-2)+Group(356,1,-2)+Group(357,1,-2)+Group(358,1,-2)+Group(359,1,-2)+Group(36,1,1)+Group(360,1,-2)+Group(361,1,-2)+Group(362,1,-2)+Group(363,1,-2)+Group(364,1,-1)+Group(365,1,-1)+Group(366,1,-1)+Group(367,1,-1)+Group(37,1,1)+Group(372,1,-2)+Group(373,1,-2)+Group(374,1,-2)+Group(375,1,-2)+Group(38,1,1)+Group(380,1,-2)+Group(381,1,-2)+Group(386,1,-1)+Group(387,1,-1)+Group(388,1,-1)+Group(389,1,-1)+Group(39,1,1)+Group(398,1,-2)+Group(399,1,-2)+Group(4,1,-1)+Group(40,1,2)+Group(400,1,-2)+Group(401,1,-2)+Group(402,1,-2)+Group(403,1,-2)+Group(404,1,-2)+Group(405,1,-2)+Group(41,1,2)+Group(410,1,-1)+Group(411,1,-1)+Group(412,1,-1)+Group(413,1,-1)+Group(414,1,-2)+Group(415,1,-2)+Group(416,1,-2)+Group(417,1,-2)+Group(418,1,-2)+Group(419,1,-2)+Group(42,1,2)+Group(420,1,-2)+Group(421,1,-2)+Group(426,1,-1)+Group(427,1,-1)+Group(428,1,-1)+Group(429,1,-1)+Group(43,1,2)+Group(430,1,-2)+Group(431,1,-2)+Group(432,1,-2)+Group(433,1,-2)+Group(434,1,-2)+Group(435,1,-2)+Group(436,1,-1)+Group(437,1,-1)+Group(438,1,-1)+Group(439,1,-1)+Group(440,1,-2)+Group(441,1,-2)+Group(442,1,-2)+Group(443,1,-2)+Group(444,1,-1)+Group(445,1,-1)+Group(446,1,-1)+Group(447,1,-1)+Group(448,1,-1/2)+Group(449,1,-1/2)+Group(450,1,-1/2)+Group(451,1,-1/2)+Group(452,1,-1)+Group(453,1,-1)+Group(454,1,-1)+Group(455,1,-1)+Group(456,1,-1/2)+Group(457,1,-1/2)+Group(458,1,-1/2)+Group(459,1,-1/2)+Group(464,1,1)+Group(465,1,1)+Group(466,1,1)+Group(467,1,1)+Group(468,1,1)+Group(469,1,1)+Group(470,1,1)+Group(471,1,1)+Group(48,1,-1)+Group(480,1,-2)+Group(481,1,-2)+Group(482,1,-2)+Group(483,1,-2)+Group(484,1,-2)+Group(485,1,-2)+Group(486,1,-1)+Group(487,1,-1)+Group(488,1,-1)+Group(489,1,-1)+Group(49,1,-1)+Group(490,1,-1)+Group(491,1,-1)+Group(492,1,-1)+Group(493,1,-1)+Group(494,1,-2)+Group(495,1,-2)+Group(496,1,-2)+Group(497,1,-2)+Group(498,1,-2)+Group(499,1,-2)+Group(5,1,1)+Group(50,1,-1)+Group(500,1,-1)+Group(501,1,-1)+Group(502,1,-1)+Group(503,1,-1)+Group(504,1,1)+Group(505,1,1)+Group(506,1,-1)+Group(507,1,-1)+Group(508,1,-1)+Group(509,1,-1)+Group(51,1,-1)+Group(512,1,-2)+Group(513,1,-2)+Group(514,1,-2)+Group(515,1,-2)+Group(516,1,-1)+Group(517,1,-1)+Group(518,1,-1)+Group(519,1,-1)+Group(522,1,-1)+Group(523,1,-1)+Group(524,1,-1)+Group(525,1,-1)+Group(526,1,-1/2)+Group(527,1,-1/2)+Group(528,1,-1/2)+Group(529,1,-1/2)+Group(534,1,-1/2)+Group(535,1,-1/2)+Group(536,1,-1/2)+Group(537,1,-1/2)+Group(538,1,-1)+Group(539,1,-1)+Group(540,1,-1)+Group(541,1,-1)+Group(542,1,1)+Group(543,1,1)+Group(552,1,-2)+Group(553,1,-2)+Group(554,1,-2)+Group(555,1,-2)+Group(556,1,-1)+Group(557,1,-1)+Group(558,1,-1)+Group(559,1,-1)+Group(56,1,-2)+Group(564,1,-1)+Group(565,1,-1)+Group(566,1,-1)+Group(567,1,-1)+Group(568,1,-2)+Group(569,1,-2)+Group(57,1,-2)+Group(570,1,-2)+Group(571,1,-2)+Group(58,1,-2)+Group(586,1,1)+Group(587,1,1)+Group(588,1,1)+Group(589,1,1)+Group(59,1,-2)+Group(590,1,2)+Group(591,1,2)+Group(592,1,2)+Group(593,1,2)+Group(594,1,1)+Group(595,1,1)+Group(596,1,1)+Group(597,1,1)+Group(598,1,2)+Group(599,1,2)+Group(6,1,1)+Group(60,1,-2)+Group(600,1,2)+Group(601,1,2)+Group(602,1,1)+Group(603,1,1)+Group(604,1,1)+Group(605,1,1)+Group(61,1,-2)+Group(610,1,1)+Group(611,1,1)+Group(612,1,1)+Group(613,1,1)+Group(614,1,2)+Group(615,1,2)+Group(616,1,2)+Group(617,1,2)+Group(62,1,-2)+Group(63,1,-2)+Group(634,1,1)+Group(635,1,1)+Group(636,1,1)+Group(637,1,1)+Group(638,1,1)+Group(639,1,1)+Group(64,1,-2)+Group(640,1,1)+Group(641,1,1)+Group(65,1,-2)+Group(658,1,1)+Group(659,1,1)+Group(66,1,-2)+Group(660,1,1)+Group(661,1,1)+Group(662,1,1)+Group(663,1,1)+Group(664,1,1)+Group(665,1,1)+Group(666,1,1)+Group(667,1,1)+Group(668,1,1)+Group(669,1,1)+Group(67,1,-2)+Group(68,1,-2)+Group(69,1,-2)+Group(694,1,1)+Group(695,1,1)+Group(696,1,1)+Group(697,1,1)+Group(7,1,-1)+Group(70,1,-2)+Group(71,1,-2)+Group(72,1,-2)+Group(73,1,-2)+Group(730,1,1)+Group(731,1,1)+Group(74,1,-2)+Group(75,1,-2)+Group(76,1,-1)+Group(77,1,-1)+Group(78,1,-1)+Group(79,1,-1)+Group(8,1,-2)+Group(9,1,2) = -426");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --symmetrize-initial-states true --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | Group(0,1,-1)+Group(1,1,-1)+Group(2,1,-1)+Group(3,1,-1) = -4");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --symmetrize-initial-states true --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"10 | Group(0,1,-2)+Group(1,1,-2)+Group(10,1,-1)+Group(11,1,-1)+Group(16,1,-1)+Group(17,1,-1)+Group(18,1,-1)+Group(19,1,-1)+Group(2,1,-2)+Group(20,1,-4)+Group(21,1,-4)+Group(22,1,-4)+Group(23,1,-4)+Group(3,1,-2)+Group(4,1,-2)+Group(5,1,-2)+Group(6,1,-2)+Group(7,1,-2)+Group(8,1,-1)+Group(9,1,-1) = -40");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g ghg QED^2==4 [{{3}} QCD=2] --symmetrize-initial-states true --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"128 | Group(0,1,-2)+Group(1,1,2)+Group(10,1,2)+Group(100,1,-2)+Group(101,1,-2)+Group(102,1,-4)+Group(103,1,-4)+Group(104,1,-4)+Group(105,1,-4)+Group(106,1,-4)+Group(107,1,-4)+Group(108,1,-2)+Group(109,1,-2)+Group(11,1,-2)+Group(110,1,-2)+Group(111,1,-2)+Group(112,1,-2)+Group(113,1,-2)+Group(114,1,-2)+Group(115,1,-2)+Group(116,1,-2)+Group(117,1,-2)+Group(118,1,-2)+Group(119,1,-2)+Group(12,1,-2)+Group(120,1,-2)+Group(121,1,-2)+Group(122,1,-2)+Group(123,1,-2)+Group(124,1,-1)+Group(125,1,-1)+Group(126,1,-1)+Group(127,1,-1)+Group(128,1,-1)+Group(129,1,-1)+Group(13,1,2)+Group(130,1,-1)+Group(131,1,-1)+Group(132,1,1)+Group(133,1,1)+Group(134,1,1)+Group(135,1,1)+Group(136,1,1)+Group(137,1,1)+Group(138,1,1)+Group(139,1,1)+Group(14,1,2)+Group(140,1,-1)+Group(141,1,-1)+Group(142,1,-1)+Group(143,1,-1)+Group(144,1,-1)+Group(145,1,-1)+Group(146,1,-1)+Group(147,1,-1)+Group(15,1,-2)+Group(152,1,-4)+Group(153,1,-4)+Group(154,1,-4)+Group(155,1,-4)+Group(16,1,1)+Group(160,1,-4)+Group(161,1,-4)+Group(162,1,-4)+Group(163,1,-4)+Group(164,1,-2)+Group(165,1,-2)+Group(166,1,-2)+Group(167,1,-2)+Group(168,1,-4)+Group(169,1,-4)+Group(17,1,1)+Group(170,1,-4)+Group(171,1,-4)+Group(172,1,2)+Group(173,1,2)+Group(174,1,2)+Group(175,1,2)+Group(176,1,-4)+Group(177,1,-4)+Group(18,1,1)+Group(182,1,-4)+Group(183,1,-4)+Group(184,1,-2)+Group(185,1,-2)+Group(186,1,-4)+Group(187,1,-4)+Group(188,1,-4)+Group(189,1,-4)+Group(19,1,1)+Group(190,1,-2)+Group(191,1,-2)+Group(192,1,-2)+Group(193,1,-2)+Group(194,1,-2)+Group(195,1,-2)+Group(196,1,-4)+Group(197,1,-4)+Group(2,1,2)+Group(202,1,-4)+Group(203,1,-4)+Group(208,1,-2)+Group(209,1,-2)+Group(210,1,-2)+Group(211,1,-2)+Group(216,1,-2)+Group(217,1,-2)+Group(218,1,-2)+Group(219,1,-2)+Group(220,1,-4)+Group(221,1,-4)+Group(222,1,-2)+Group(223,1,-2)+Group(228,1,-2)+Group(229,1,-2)+Group(230,1,-2)+Group(231,1,-2)+Group(232,1,-4)+Group(233,1,-4)+Group(234,1,-4)+Group(235,1,-4)+Group(236,1,-4)+Group(237,1,-4)+Group(238,1,-2)+Group(239,1,-2)+Group(24,1,1)+Group(240,1,-2)+Group(241,1,-2)+Group(242,1,-2)+Group(243,1,-2)+Group(244,1,-2)+Group(245,1,-2)+Group(246,1,-1)+Group(247,1,-1)+Group(248,1,-1)+Group(249,1,-1)+Group(25,1,1)+Group(250,1,-1)+Group(251,1,-1)+Group(252,1,-1)+Group(253,1,-1)+Group(254,1,-1)+Group(255,1,-1)+Group(256,1,-1)+Group(257,1,-1)+Group(258,1,-1)+Group(259,1,-1)+Group(26,1,1)+Group(260,1,-1)+Group(261,1,-1)+Group(262,1,-4)+Group(263,1,-4)+Group(264,1,-4)+Group(265,1,-4)+Group(266,1,-2)+Group(267,1,-2)+Group(268,1,-1)+Group(269,1,-1)+Group(27,1,1)+Group(270,1,-1)+Group(271,1,-1)+Group(272,1,1)+Group(273,1,1)+Group(274,1,-2)+Group(275,1,-2)+Group(276,1,-2)+Group(277,1,-2)+Group(278,1,-2)+Group(279,1,-2)+Group(28,1,4)+Group(280,1,-2)+Group(281,1,-2)+Group(286,1,-1/2)+Group(287,1,-1/2)+Group(288,1,-1/2)+Group(289,1,-1/2)+Group(29,1,4)+Group(294,1,-1/2)+Group(295,1,-1/2)+Group(296,1,-1/2)+Group(297,1,-1/2)+Group(298,1,-2)+Group(299,1,-2)+Group(3,1,-2)+Group(30,1,4)+Group(300,1,-2)+Group(301,1,-2)+Group(302,1,-1)+Group(303,1,-1)+Group(304,1,-1)+Group(305,1,-1)+Group(31,1,4)+Group(310,1,-1)+Group(311,1,-1)+Group(312,1,-1)+Group(313,1,-1)+Group(314,1,-4)+Group(315,1,-4)+Group(316,1,-4)+Group(317,1,-4)+Group(32,1,-2)+Group(33,1,-2)+Group(34,1,-4)+Group(35,1,-4)+Group(36,1,-4)+Group(37,1,-4)+Group(38,1,-2)+Group(39,1,-2)+Group(4,1,-2)+Group(40,1,-4)+Group(41,1,-4)+Group(42,1,-4)+Group(43,1,-4)+Group(44,1,-2)+Group(45,1,-2)+Group(46,1,-2)+Group(47,1,-2)+Group(5,1,2)+Group(52,1,-4)+Group(53,1,-4)+Group(54,1,-4)+Group(55,1,-4)+Group(56,1,-4)+Group(57,1,-4)+Group(58,1,-4)+Group(59,1,-4)+Group(6,1,2)+Group(60,1,-4)+Group(61,1,-4)+Group(66,1,-4)+Group(67,1,-4)+Group(68,1,-4)+Group(69,1,-4)+Group(7,1,-2)+Group(70,1,-4)+Group(71,1,-4)+Group(72,1,-4)+Group(73,1,-4)+Group(74,1,-4)+Group(75,1,-4)+Group(76,1,-4)+Group(77,1,-4)+Group(78,1,-4)+Group(79,1,-4)+Group(8,1,-2)+Group(80,1,-2)+Group(81,1,-2)+Group(82,1,-2)+Group(83,1,-2)+Group(84,1,-2)+Group(85,1,-2)+Group(86,1,-2)+Group(87,1,-2)+Group(88,1,-2)+Group(89,1,-2)+Group(9,1,2)+Group(90,1,-4)+Group(91,1,-4)+Group(92,1,-2)+Group(93,1,-2)+Group(98,1,-2)+Group(99,1,-2) = -502");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_aa_ttx_cross_section_slow_filter() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --numerator-grouping only_detect_zeroes --max-multiplicity-for-fast-cut-filter 0",false)?,@"4 | -4 = -4");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --numerator-grouping only_detect_zeroes --max-multiplicity-for-fast-cut-filter 0",false)?,@"40 | -40 = -40");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] --numerator-grouping only_detect_zeroes --max-multiplicity-for-fast-cut-filter 0",false)?,@"874 | -266 = -266");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_vacuum_amplitude_kaapo() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm-default.json"))?;

    // 4-loop vaccuum contribution to the neutron start equation of state
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g d d~ ghG ghG~ [{4}] --numerator-grouping only_detect_zeroes --number-of-factorized-loop-subtopologies 1 1000 --number-of-fermion-loops 1 1000 --filter-snails false --filter-selfenergies false --filter-tadpoles false --max-n-bridges 0",false)?,@"52 | -44/3 = -44/3");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_vacuum_amplitude_generation() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Including all graphs
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{2}] --numerator-grouping only_detect_zeroes --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 -1 ",false)?,@"6 | -19/24 = -19/24");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 -1",false)?,@"5 | 5/24+Group(0,1,-1/2)+Group(5,1,-1/2) = -19/24");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{3}] --numerator-grouping only_detect_zeroes --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 -1",false)?,@"36 | -83/48 = -83/48");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{3}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 -1",false)?,@"27 | 7/16+Group(0,1,-1/4)+Group(1,1,-1/4)+Group(11,1,1/2)+Group(16,1,-1/2)+Group(19,1,-1/2)+Group(2,1,-1/4)+Group(3,1,-1/3)+Group(30,1,-1/4)+Group(31,1,-1/4)+Group(32,1,-1/3)+Group(33,1,-1/4)+Group(34,1,1/4)+Group(35,1,-1/2)+Group(4,1,1/4)+Group(5,1,-1/2)+Group(6,1,1/2)+Group(7,1,1/2) = -83/48");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{4}] --numerator-grouping only_detect_zeroes --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 -1",false)?,@"264 | -143/24 = -143/24");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{4}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 -1",false)?,@"171 | 13/8+Group(0,1,-1/8)+Group(1,1,-1/2)+Group(10,1,-1/8)+Group(100,1,1)+Group(101,1,1)+Group(105,1,1)+Group(107,1,1/2)+Group(11,1,-1/3)+Group(111,1,1/2)+Group(115,1,-1/2)+Group(119,1,-1/2)+Group(12,1,-1/2)+Group(126,1,1/6)+Group(127,7/2,1/6)+Group(13,1,-1/2)+Group(134,1,-1/2)+Group(135,1,-1/2)+Group(136,1,-1)+Group(137,1,-1/2)+Group(138,1,-1/4)+Group(139,1,-1/2)+Group(14,1,-1/2)+Group(140,1,-1)+Group(141,1,1/2)+Group(15,1,-1/2)+Group(153,1,1)+Group(154,1,-1/2)+Group(155,1,-1/2)+Group(156,1,-1/2)+Group(157,1,-1/4)+Group(158,1,-1)+Group(159,1,-1)+Group(16,1,-1/2)+Group(160,1,-1/2)+Group(161,1,1/2)+Group(162,1,-1/3)+Group(163,-1,-1/3)+Group(169,1,-1/3)+Group(17,1,-1/2)+Group(170,1,-1/3)+Group(172,1,-1)+Group(175,1,-1)+Group(177,1,-1/2)+Group(18,1,-1/4)+Group(180,1,-1/2)+Group(184,1,-1/2)+Group(187,1,-1/2)+Group(19,1,-1/4)+Group(2,1,-1/4)+Group(21,1,-1/4)+Group(22,1,-1)+Group(228,1,1)+Group(229,1,1/2)+Group(23,1,1/4)+Group(230,1,1)+Group(231,1,-1/8)+Group(232,1,-1/2)+Group(233,1,-1/4)+Group(234,1,-1/8)+Group(235,1,-1/4)+Group(236,1,-1/8)+Group(237,1,-1/4)+Group(238,1,-1/8)+Group(239,1,-1/12)+Group(24,1,1/4)+Group(240,1,-1/4)+Group(241,1,-1/8)+Group(242,1,-1/3)+Group(243,1,-1/2)+Group(244,1,-1/2)+Group(245,1,-1/2)+Group(246,1,-1/2)+Group(247,1,-1/2)+Group(248,1,-1/2)+Group(249,1,-1/4)+Group(25,1,1/4)+Group(250,1,-1/4)+Group(252,1,-1/4)+Group(253,1,-1)+Group(254,1,1/4)+Group(255,1,1/4)+Group(256,1,1/4)+Group(257,1,1/8)+Group(258,1,-1)+Group(259,1,-1/6)+Group(26,1,1/8)+Group(260,1,1)+Group(261,1,-1/2)+Group(262,1,-1)+Group(263,1,1/6)+Group(264,7/2,1/6)+Group(265,1,1/2)+Group(266,1,1)+Group(267,1,-1/3)+Group(268,1,-1/2)+Group(269,1,-1/6)+Group(27,1,-1)+Group(28,1,-1/6)+Group(29,1,1)+Group(3,1,-1/8)+Group(30,1,-1/2)+Group(31,1,-1/2)+Group(32,1,1/6)+Group(33,7/2,1/6)+Group(34,1,-1)+Group(35,1,1/2)+Group(36,1,1)+Group(37,1,-1/3)+Group(38,1,-1/2)+Group(39,1,-1/6)+Group(4,1,-1/4)+Group(40,1,1/2)+Group(41,1,1/2)+Group(42,1,1/2)+Group(43,1,1/4)+Group(44,1,1)+Group(45,1,-1/2)+Group(46,7/2,1/3)+Group(47,1,1/3)+Group(48,1,1/2)+Group(49,1,1)+Group(5,1,-1/8)+Group(50,1,1/2)+Group(51,1,1/2)+Group(52,1,1)+Group(53,1,1/2)+Group(54,1,1/4)+Group(55,1,1/2)+Group(56,1,1)+Group(57,1,-1/2)+Group(58,1,1)+Group(6,1,-1/4)+Group(61,1,1)+Group(7,1,-1/8)+Group(79,1,-1)+Group(8,1,-1/12)+Group(80,1,1/2)+Group(81,1,1/2)+Group(82,1,1/2)+Group(83,1,1/4)+Group(84,1,1)+Group(85,1,1)+Group(86,1,1/2)+Group(87,1,-1/2)+Group(88,1,1/3)+Group(89,2/7,1/3)+Group(9,1,-1/4)+Group(90,1,1)+Group(98,1,1/3)+Group(99,2/7,1/3) = -311/168");

    // // Only non-factorizable ones
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{2}] --numerator-grouping only_detect_zeroes --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1 1",false)?,@"5 | -11/12 = -11/12");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1 1",false)?,@"4 | 1/12+Group(0,1,-1/2)+Group(4,1,-1/2) = -11/12");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{3}] --numerator-grouping only_detect_zeroes --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1 1",false)?,@"30 | -17/12 = -17/12");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{3}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1 1",false)?,@"22 | 1/4+Group(0,1,-1/4)+Group(1,1,-1/4)+Group(14,1,-1/2)+Group(16,1,-1/2)+Group(2,1,-1/3)+Group(25,1,-1/4)+Group(26,1,-1/3)+Group(27,1,-1/4)+Group(28,1,1/4)+Group(29,1,-1/2)+Group(3,1,1/4)+Group(4,1,-1/2)+Group(5,1,1/2)+Group(6,1,1/2)+Group(9,1,1/2) = -17/12");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{3}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 1 --number-of-factorized-loop-subtopologies 1 1",true)?,@"22 | 1/4+Group(0,1,-1/4)+Group(1,1,-1/4)+Group(14,1,-1/2)+Group(16,1,-1/2)+Group(2,1,-1/3)+Group(25,1,-1/4)+Group(26,1,-1/3)+Group(27,1,-1/4)+Group(28,1,1/4)+Group(29,1,-1/2)+Group(3,1,1/4)+Group(4,1,-1/2)+Group(5,1,1/2)+Group(6,1,1/2)+Group(9,1,1/2) = -17/12");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{4}] --numerator-grouping only_detect_zeroes --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1 1",false)?,@"200 | -199/48 = -199/48");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{4}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1 1",false)?,@"123 | 11/16+Group(0,1,-1/8)+Group(1,1,-1/2)+Group(10,1,-1/2)+Group(106,1,-1/2)+Group(107,1,-1)+Group(108,1,-1/2)+Group(109,1,-1/2)+Group(11,1,-1/4)+Group(110,1,-1)+Group(111,1,1/2)+Group(118,1,1)+Group(119,1,-1/2)+Group(120,1,-1/2)+Group(121,1,-1)+Group(122,1,-1)+Group(123,1,-1/2)+Group(124,1,1/2)+Group(125,1,-1/3)+Group(126,-1,-1/3)+Group(13,1,-1/4)+Group(131,1,-1/3)+Group(132,1,-1/3)+Group(134,1,-1)+Group(136,1,-1)+Group(138,1,-1/2)+Group(14,1,-1)+Group(140,1,-1/2)+Group(144,1,-1/2)+Group(146,1,-1/2)+Group(15,1,1/4)+Group(16,1,1/4)+Group(17,1,-1)+Group(174,1,1)+Group(175,1,1/2)+Group(176,1,1)+Group(177,1,-1/8)+Group(178,1,-1/2)+Group(179,1,-1/4)+Group(18,1,-1/6)+Group(180,1,-1/4)+Group(181,1,-1/8)+Group(182,1,-1/12)+Group(183,1,-1/3)+Group(184,1,-1/2)+Group(185,1,-1/2)+Group(186,1,-1/2)+Group(187,1,-1/2)+Group(188,1,-1/4)+Group(19,1,1)+Group(190,1,-1/4)+Group(191,1,-1)+Group(192,1,1/4)+Group(193,1,1/4)+Group(194,1,-1)+Group(195,1,-1/6)+Group(196,1,1)+Group(197,1,-1/2)+Group(198,1,-1)+Group(199,1,1/6)+Group(2,1,-1/4)+Group(20,1,-1/2)+Group(200,7/2,1/6)+Group(201,1,1/2)+Group(202,1,1)+Group(203,1,-1/3)+Group(204,1,-1/2)+Group(205,1,-1/6)+Group(21,1,-1/2)+Group(22,1,1/6)+Group(23,7/2,1/6)+Group(24,1,-1)+Group(25,1,1/2)+Group(26,1,1)+Group(27,1,-1/3)+Group(28,1,-1/2)+Group(29,1,-1/6)+Group(3,1,-1/4)+Group(30,1,1/2)+Group(31,1,1/2)+Group(32,1,1)+Group(33,1,-1/2)+Group(34,7/2,1/3)+Group(35,1,1/3)+Group(36,1,1/2)+Group(37,1,1)+Group(38,1,1/2)+Group(39,1,1)+Group(4,1,-1/8)+Group(40,1,1/2)+Group(41,1,1/2)+Group(42,1,1)+Group(43,1,-1/2)+Group(44,1,1)+Group(46,1,1)+Group(5,1,-1/12)+Group(58,1,-1)+Group(59,1,1/2)+Group(6,1,-1/3)+Group(60,1,1/2)+Group(61,1,1)+Group(62,1,1)+Group(63,1,1/2)+Group(64,1,-1/2)+Group(65,1,1/3)+Group(66,2/7,1/3)+Group(67,1,1)+Group(7,1,-1/2)+Group(74,1,1/3)+Group(75,2/7,1/3)+Group(76,1,1)+Group(77,1,1)+Group(8,1,-1/2)+Group(80,1,1)+Group(82,1,1/2)+Group(85,1,1/2)+Group(89,1,-1/2)+Group(9,1,-1/2)+Group(92,1,-1/2)+Group(98,1,1/6)+Group(99,7/2,1/6) = -13/336");

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_slow_vacuum_amplitude_generation() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // CURRENTLY BUGGED BECAUSE OF SYMBOLICA INCORRECT SYMMETRY FACTORS. POSSIBLY OTHER SCENARIOS BUGGED TOO.
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{5}] --numerator-grouping only_detect_zeroes --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1 1",false)?,@"2560 | -5785/384");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{5}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1 1",false)?,@"1440 | 233015/51072");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_dis_isr() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // 2>N amplitude, no symmetrization
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d > e- d | e- a d g ghg QED==2 [{0}] --numerator-grouping only_detect_zeroes",false)?,@"1 | -1");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d > e- d | e- a d g ghg QED==2 [{1}] --numerator-grouping only_detect_zeroes",false)?,@"1 | -1");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d > e- d | e- a d g ghg QED==2 [{2}] --numerator-grouping only_detect_zeroes",false)?,@"13 | -13/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d > e- d | e- a d g ghg QED==2 [{3}] --numerator-grouping only_detect_zeroes",false)?,@"241 | -194/3");

    // 2>N amplitude, with symmetrization
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d > e- d | e- a d g ghg QED==2 [{0}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"1 | -1");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d > e- d | e- a d g ghg QED==2 [{1}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"1 | -1");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d > e- d | e- a d g ghg QED==2 [{2}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"11 | -17/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d > e- d | e- a d g ghg QED==2 [{3}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"191 | -2183/21");

    // 3>N no symmetrization
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d g > e- d | e- a d g ghg QED==2 [{0}] --numerator-grouping only_detect_zeroes",false)?,@"2 | -2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d g > e- d | e- a d g ghg QED==2 [{1}] --numerator-grouping only_detect_zeroes",false)?,@"13 | -9");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d g > e- d | e- a d g ghg QED==2 [{2}] --numerator-grouping only_detect_zeroes",false)?,@"225 | -169/2");

    // 3>N with symmetrization
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d g > e- d | e- a d g ghg QED==2 [{0}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d g > e- d | e- a d g ghg QED==2 [{1}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"11 | -11");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d g > e- d | e- a d g ghg QED==2 [{2}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"181 | -842/7");

    // 4>N no symmetrization
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d u u~ > e- d | e- a d u g ghg QED==2 [{0}] --numerator-grouping only_detect_zeroes",false)?,@"4 | 4");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d u u~ > e- d | e- a d u g ghg QED==2 [{1}] --numerator-grouping only_detect_zeroes",false)?,@"62 | 36");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d u u~ > e- d | e- a d u g ghg QED==2 [{2}] --numerator-grouping only_detect_zeroes",false)?,@"1428 | 1298/3");

    // 4>N with symmetrization
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d u u~ > e- d | e- a d u g ghg QED==2 [{0}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"4 | 4");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d u u~ > e- d | e- a d u g ghg QED==2 [{1}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"54 | 40");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e- d u u~ > e- d | e- a d u g ghg QED==2 [{2}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"1087 | 22357/42");

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_slow_generate_a_qqh() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "amp", "a > b b~ h | b h a ghg g QED==4 [{5} QCD=3] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2763 | -38712/7");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_amplitude_1l_sm_jets() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Targets confirmed by MadGraph
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"20 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"192 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > g g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"384 | -186");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "g g > g g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"1139 | 171/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"2284 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"10 | -11/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"96 | -119/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"1142 | -698");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ d d~ | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"414 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ d d~ | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"1242 | 0");

    Ok(())
}

// Currently bugged with spenso network parsing for "--numerator-grouping group_identical_graphs_up_to_sign"
#[test]
#[rustfmt::skip]
fn test_generate_amplitude_1l_sm_jets_with_grouping() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Targets confirmed by MadGraph
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"18 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"9 | -11/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"88 | -119/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"176 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > g g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"341 | -186");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "g g > g g g | g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"671 | 1107/2");
    // For fun, symmetrize it all
    assert_snapshot!(feyngen_str(&mut cli, "amp", "g g > g g g | g ghg a QED==0 [QCD=1] --symmetrize-initial-states --symmetrize-final-states --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"109 | 1107/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "g g > g g g | g ghg a QED==0 [QCD=1] --symmetrize-left-right-states true --symmetrize-initial-states --symmetrize-final-states --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"25 | 1107/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "g g > g g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"905 | 171/2");

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_slow_generate_amplitude_1l_sm_jets() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Targets confirmed by MadGraph
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ d d~ g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"5424 | 0");

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "veryslow"]
fn test_very_slow_generate_amplitude_1l_sm_jets() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Targets confirmed by MadGraph
    assert_snapshot!(feyngen_str(&mut cli, "amp", "g g > g g g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"14875 | 1375");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ g g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"32074 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ g g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"16037 | -18993/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ d d~ g | u d g ghg a QED==0 [QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"16272 | 0");

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "veryslow"]
fn test_slow_generate_amplitude_1l_sm_jets_with_grouping() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Targets confirmed by MadGraph
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"2090 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"1045 | -698");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ d d~ | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"380 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ d d~ | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"1140 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "g g > g g g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"11850 | 1375");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ g g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"29210 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ g g g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"14605 | -18993/2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > d d~ d d~ g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"15030 | 0");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "d d~ > u u~ d d~ g | u d g ghg a QED==0 [QCD=1] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"5010 | 0");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_vaccuum_amplitude_generation_full_sm() -> Result<()> {
    let mut cli = get_test_cli(None,get_tests_workspace_path().join("feyn_gen_generation_test"),Some("feyngen".to_string()),true)?;

    // Choose the model to consider
    cli.run_command(&format!("import model sm-full.json"))?;

    // Including all graphs
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{2}] --numerator-grouping only_detect_zeroes --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1",false)?,@"448 | TBD");

    // Compare various numerator comparison options below:
    // a) Using only symbolic tensor canonization of the numerator and without substituting color group factors, some groupings are missed (i.e. 293 graphs instead of 243 when grouping as much as possible)
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 --compare-canonized-numerator true --number-of-samples-for-numerator-comparisons 0 --fully-numerical-substitution-when-comparing-numerators",false)?,@"281 | (UFO::complexconjugate(UFO::CKM1x1))^(-1)*-1*UFO::CKM1x1^(-1)*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM1x2)+(UFO::complexconjugate(UFO::CKM2x1))^(-1)*-1*UFO::CKM2x1^(-1)*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM2x2)+(UFO::complexconjugate(UFO::CKM3x1))^(-1)*-1*UFO::CKM3x1^(-1)*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM3x2)+-1*UFO::I2x12^(-1)*UFO::I2x22*UFO::I3x21^(-1)*UFO::I3x22+-1*UFO::I2x13^(-1)*UFO::I2x23*UFO::I3x31^(-1)*UFO::I3x32+-27/2*UFO::G^2*UFO::ee^(-2)+319/24");
    // b) Using only symbolic tensor canonization of the numerator but this time substituting color group factors allows for more groupings, but still not the maximum (i.e. 281 graphs instead of 243 when grouping as much as possible)


    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 --compare-canonized-numerator true --number-of-samples-for-numerator-comparisons 0 --fully-numerical-substitution-when-comparing-numerators false",false)?,@"281 | -UFO::CKM1x1^-1*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22-UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32-81/8*UFO::G^2*UFO::ee^-2*UFO::Nc*UFO::TR+81/8*UFO::G^2*UFO::ee^-2*UFO::Nc^-1*UFO::TR+319/24");

    // c) Now instead rely on numerical sample evaluations but still with symbolic color group factors
    // Generation below yields "Could not numerically evaluate numerator: more than one node in the graph" with symbolic EW parameters
    // assert_snapshot!(generate_graphs_and_count_to_str(&mut cli, "amp", "{} > {} [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 --compare-canonized-numerator false --number-of-samples-for-numerator-comparisons 3 --no-fully-numerical-substitution-when-comparing-numerators",false)?,@"? | ?");

    // d) Fully numerical evaluations of the numerator, yielding the largest amount of groupings (i.e. 243 graphs)
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 --compare-canonized-numerator false --number-of-samples-for-numerator-comparisons 3 --fully-numerical-substitution-when-comparing-numerators",false)?,@"243 | -27/2*UFO::G^2*UFO::ee^-2-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22-UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32+313/24");

    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1",false)?,@"243 | -27/2*UFO::G^2*UFO::ee^-2-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22-UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32+313/24");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{3}] --numerator-grouping only_detect_zeroes --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1",false)?,@"36362 | -68");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{3}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1",false)?,@"13869 | 81*UFO::G^2*UFO::ee^-2-1053/16*UFO::G^4*UFO::ee^-4-31/2*UFO::cw^2*(-UFO::cw^2+UFO::sw^2)^-1-1/2*UFO::cw^4*(-1/2*UFO::cw^4-1/2*UFO::sw^4+UFO::cw^2*UFO::sw^2)^-1-31/2*UFO::sw^2*(-UFO::cw^2+UFO::sw^2)^-1-1/2*UFO::sw^4*(-1/2*UFO::cw^4-1/2*UFO::sw^4+UFO::cw^2*UFO::sw^2)^-1-UFO::cw^2*UFO::sw^2*(-1/2*UFO::cw^4-1/2*UFO::sw^4+UFO::cw^2*UFO::sw^2)^-1-5/2*UFO::CKM2x1*UFO::I2x12*UFO::I3x21^-1*UFO::complexconjugate(UFO::CKM2x1)^-1-5/2*UFO::CKM2x2*UFO::I2x22*UFO::I3x21^-1*UFO::complexconjugate(UFO::CKM2x1)^-1+1/4*UFO::CKM2x2*UFO::I2x22*UFO::cw^2*(-UFO::CKM2x1*UFO::I2x12*UFO::cw^2+UFO::CKM2x1*UFO::I2x12*UFO::sw^2)^-1+1/4*UFO::CKM2x2*UFO::I2x22*UFO::sw^2*(-UFO::CKM2x1*UFO::I2x12*UFO::cw^2+UFO::CKM2x1*UFO::I2x12*UFO::sw^2)^-1-5/2*UFO::CKM3x1*UFO::I2x13*UFO::I3x31^-1*UFO::complexconjugate(UFO::CKM3x1)^-1-5/2*UFO::CKM3x2*UFO::I2x23*UFO::I3x31^-1*UFO::complexconjugate(UFO::CKM3x1)^-1+1/4*UFO::CKM3x2*UFO::I2x23*UFO::cw^2*(-UFO::CKM3x1*UFO::I2x13*UFO::cw^2+UFO::CKM3x1*UFO::I2x13*UFO::sw^2)^-1+1/4*UFO::CKM3x2*UFO::I2x23*UFO::sw^2*(-UFO::CKM3x1*UFO::I2x13*UFO::cw^2+UFO::CKM3x1*UFO::I2x13*UFO::sw^2)^-1+1/4*UFO::I3x22*UFO::cw^2*(-UFO::I3x21*UFO::cw^2*UFO::complexconjugate(UFO::CKM2x1)+UFO::I3x21*UFO::sw^2*UFO::complexconjugate(UFO::CKM2x1))^-1*UFO::complexconjugate(UFO::CKM2x2)+1/4*UFO::I3x22*UFO::sw^2*(-UFO::I3x21*UFO::cw^2*UFO::complexconjugate(UFO::CKM2x1)+UFO::I3x21*UFO::sw^2*UFO::complexconjugate(UFO::CKM2x1))^-1*UFO::complexconjugate(UFO::CKM2x2)+1/4*UFO::I3x32*UFO::cw^2*(-UFO::I3x31*UFO::cw^2*UFO::complexconjugate(UFO::CKM3x1)+UFO::I3x31*UFO::sw^2*UFO::complexconjugate(UFO::CKM3x1))^-1*UFO::complexconjugate(UFO::CKM3x2)+1/4*UFO::I3x32*UFO::sw^2*(-UFO::I3x31*UFO::cw^2*UFO::complexconjugate(UFO::CKM3x1)+UFO::I3x31*UFO::sw^2*UFO::complexconjugate(UFO::CKM3x1))^-1*UFO::complexconjugate(UFO::CKM3x2)-5/2*UFO::CKM1x3*UFO::I1x31*UFO::I4x13^-1*UFO::complexconjugate(UFO::CKM1x3)^-1-1/2*UFO::CKM1x1^-2*UFO::CKM1x2^2*UFO::complexconjugate(UFO::CKM1x1)^-2*UFO::complexconjugate(UFO::CKM1x2)^2-2*UFO::CKM1x1^-1*UFO::CKM1x2*UFO::I2x12^-1*UFO::I2x22-2*UFO::CKM1x1^-1*UFO::CKM1x2*UFO::I2x13^-1*UFO::I2x23-1/2*UFO::CKM1x1^-1*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-5*UFO::CKM1x1^-1*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-5*UFO::CKM1x1^-1*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-1/2*UFO::CKM2x1^-2*UFO::CKM2x2^2*UFO::complexconjugate(UFO::CKM2x1)^-2*UFO::complexconjugate(UFO::CKM2x2)^2+31/4*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::I2x12^-1*UFO::I2x22-2*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::I2x13^-1*UFO::I2x23-5*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)+5/2*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-6*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-1/2*UFO::CKM3x1^-2*UFO::CKM3x2^2*UFO::complexconjugate(UFO::CKM3x1)^-2*UFO::complexconjugate(UFO::CKM3x2)^2-2*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::I2x12^-1*UFO::I2x22+31/4*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::I2x13^-1*UFO::I2x23-5*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-6*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)+5/2*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-1/2*UFO::I2x12^-2*UFO::I2x22^2*UFO::I3x21^-2*UFO::I3x22^2+7/2*UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22-6*UFO::I2x12^-1*UFO::I2x22*UFO::I3x31^-1*UFO::I3x32-6*UFO::I3x21^-1*UFO::I3x22*UFO::I2x13^-1*UFO::I2x23-2*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)+11/2*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-2*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-1/2*UFO::I2x13^-2*UFO::I2x23^2*UFO::I3x31^-2*UFO::I3x32^2+7/2*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32-2*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-2*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)+11/2*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-9*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::G^2*UFO::ee^-2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-9*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::G^2*UFO::ee^-2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-9*UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22*UFO::G^2*UFO::ee^-2-9*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32*UFO::G^2*UFO::ee^-2-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-UFO::CKM2x1^-1*UFO::CKM2x2*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::CKM2x1^-1*UFO::CKM2x2*UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-UFO::CKM2x1^-1*UFO::CKM2x2*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-UFO::CKM3x1^-1*UFO::CKM3x2*UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::CKM3x1^-1*UFO::CKM3x2*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32+4525/16");

    // Only non-factorizable ones
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{2}] --numerator-grouping only_detect_zeroes --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"92 | -67/3");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"69 | -27/2*UFO::G^2*UFO::ee^-2-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22-UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32-163/12");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{3}] --numerator-grouping only_detect_zeroes --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"3085 | -1687/3");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} [{3}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"2009 | -981/4*UFO::G^2*UFO::ee^-2-1053/16*UFO::G^4*UFO::ee^-4-9*UFO::cw^2*(-UFO::cw^2+UFO::sw^2)^-1-1/2*UFO::cw^4*(-1/2*UFO::cw^4-1/2*UFO::sw^4+UFO::cw^2*UFO::sw^2)^-1-9*UFO::sw^2*(-UFO::cw^2+UFO::sw^2)^-1-1/2*UFO::sw^4*(-1/2*UFO::cw^4-1/2*UFO::sw^4+UFO::cw^2*UFO::sw^2)^-1-UFO::cw^2*UFO::sw^2*(-1/2*UFO::cw^4-1/2*UFO::sw^4+UFO::cw^2*UFO::sw^2)^-1-5/2*UFO::CKM2x1*UFO::I2x12*UFO::I3x21^-1*UFO::complexconjugate(UFO::CKM2x1)^-1-5/2*UFO::CKM2x2*UFO::I2x22*UFO::I3x21^-1*UFO::complexconjugate(UFO::CKM2x1)^-1+1/4*UFO::CKM2x2*UFO::I2x22*UFO::cw^2*(-UFO::CKM2x1*UFO::I2x12*UFO::cw^2+UFO::CKM2x1*UFO::I2x12*UFO::sw^2)^-1+1/4*UFO::CKM2x2*UFO::I2x22*UFO::sw^2*(-UFO::CKM2x1*UFO::I2x12*UFO::cw^2+UFO::CKM2x1*UFO::I2x12*UFO::sw^2)^-1-5/2*UFO::CKM3x1*UFO::I2x13*UFO::I3x31^-1*UFO::complexconjugate(UFO::CKM3x1)^-1-5/2*UFO::CKM3x2*UFO::I2x23*UFO::I3x31^-1*UFO::complexconjugate(UFO::CKM3x1)^-1+1/4*UFO::CKM3x2*UFO::I2x23*UFO::cw^2*(-UFO::CKM3x1*UFO::I2x13*UFO::cw^2+UFO::CKM3x1*UFO::I2x13*UFO::sw^2)^-1+1/4*UFO::CKM3x2*UFO::I2x23*UFO::sw^2*(-UFO::CKM3x1*UFO::I2x13*UFO::cw^2+UFO::CKM3x1*UFO::I2x13*UFO::sw^2)^-1+1/4*UFO::I3x22*UFO::cw^2*(-UFO::I3x21*UFO::cw^2*UFO::complexconjugate(UFO::CKM2x1)+UFO::I3x21*UFO::sw^2*UFO::complexconjugate(UFO::CKM2x1))^-1*UFO::complexconjugate(UFO::CKM2x2)+1/4*UFO::I3x22*UFO::sw^2*(-UFO::I3x21*UFO::cw^2*UFO::complexconjugate(UFO::CKM2x1)+UFO::I3x21*UFO::sw^2*UFO::complexconjugate(UFO::CKM2x1))^-1*UFO::complexconjugate(UFO::CKM2x2)+1/4*UFO::I3x32*UFO::cw^2*(-UFO::I3x31*UFO::cw^2*UFO::complexconjugate(UFO::CKM3x1)+UFO::I3x31*UFO::sw^2*UFO::complexconjugate(UFO::CKM3x1))^-1*UFO::complexconjugate(UFO::CKM3x2)+1/4*UFO::I3x32*UFO::sw^2*(-UFO::I3x31*UFO::cw^2*UFO::complexconjugate(UFO::CKM3x1)+UFO::I3x31*UFO::sw^2*UFO::complexconjugate(UFO::CKM3x1))^-1*UFO::complexconjugate(UFO::CKM3x2)-5/2*UFO::CKM1x3*UFO::I1x31*UFO::I4x13^-1*UFO::complexconjugate(UFO::CKM1x3)^-1-1/2*UFO::CKM1x1^-2*UFO::CKM1x2^2*UFO::complexconjugate(UFO::CKM1x1)^-2*UFO::complexconjugate(UFO::CKM1x2)^2-2*UFO::CKM1x1^-1*UFO::CKM1x2*UFO::I2x12^-1*UFO::I2x22-2*UFO::CKM1x1^-1*UFO::CKM1x2*UFO::I2x13^-1*UFO::I2x23-63/4*UFO::CKM1x1^-1*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-5*UFO::CKM1x1^-1*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-5*UFO::CKM1x1^-1*UFO::CKM1x2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-1/2*UFO::CKM2x1^-2*UFO::CKM2x2^2*UFO::complexconjugate(UFO::CKM2x1)^-2*UFO::complexconjugate(UFO::CKM2x2)^2+7/2*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::I2x12^-1*UFO::I2x22-2*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::I2x13^-1*UFO::I2x23-5*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-81/4*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-6*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-1/2*UFO::CKM3x1^-2*UFO::CKM3x2^2*UFO::complexconjugate(UFO::CKM3x1)^-2*UFO::complexconjugate(UFO::CKM3x2)^2-2*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::I2x12^-1*UFO::I2x22+7/2*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::I2x13^-1*UFO::I2x23-5*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-6*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-81/4*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-1/2*UFO::I2x12^-2*UFO::I2x22^2*UFO::I3x21^-2*UFO::I3x22^2-77/4*UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22-6*UFO::I2x12^-1*UFO::I2x22*UFO::I3x31^-1*UFO::I3x32-6*UFO::I3x21^-1*UFO::I3x22*UFO::I2x13^-1*UFO::I2x23-2*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)+5/4*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-2*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-1/2*UFO::I2x13^-2*UFO::I2x23^2*UFO::I3x31^-2*UFO::I3x32^2-77/4*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32-2*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-2*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)+5/4*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-9*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::G^2*UFO::ee^-2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-9*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::G^2*UFO::ee^-2*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-9*UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22*UFO::G^2*UFO::ee^-2-9*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32*UFO::G^2*UFO::ee^-2-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::CKM2x1^-1*UFO::CKM2x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-UFO::CKM1x1^-1*UFO::CKM1x2*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM1x1)^-1*UFO::complexconjugate(UFO::CKM1x2)-UFO::CKM2x1^-1*UFO::CKM2x2*UFO::CKM3x1^-1*UFO::CKM3x2*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::CKM2x1^-1*UFO::CKM2x2*UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-UFO::CKM2x1^-1*UFO::CKM2x2*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM2x1)^-1*UFO::complexconjugate(UFO::CKM2x2)-UFO::CKM3x1^-1*UFO::CKM3x2*UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::CKM3x1^-1*UFO::CKM3x2*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32*UFO::complexconjugate(UFO::CKM3x1)^-1*UFO::complexconjugate(UFO::CKM3x2)-UFO::I2x12^-1*UFO::I2x22*UFO::I3x21^-1*UFO::I3x22*UFO::I2x13^-1*UFO::I2x23*UFO::I3x31^-1*UFO::I3x32-25787/96");

    Ok(())
}
