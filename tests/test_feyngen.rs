use bincode::de;
use color_eyre::Result;

mod test_utils;
use gammaloop_api::{commands::save, state::SyncSettings};
use gammalooprs::feyngen::diagram_generator::evaluate_overall_factor;
use gammalooprs::processes::ProcessCollection;
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
        overall_factor_sum.to_canonically_ordered_string(CanonicalOrderingSettings {
            include_namespace: false,
            include_attributes: false,
            ..Default::default()
        }),
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
      overall_factor = "(AutG(1))^(-1)*AntiFermionSpinSumSign(-1)*ExternalFermionOrderingSign(-1)*InternalFermionLoopSign(-1)";
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
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ / z QED^2==4 [{{1}} QCD] --numerator-grouping no_grouping",false)?,@"1 | (AutG(1))^(-1)*AntiFermionSpinSumSign(-1)*ExternalFermionOrderingSign(-1)*InternalFermionLoopSign(-1) = -1");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ | e- a d g QED^2==4 [{{2}} QCD] --numerator-grouping only_detect_zeroes",false)?,@"3 | (AutG(1))^(-1)*3*AntiFermionSpinSumSign(-1)*ExternalFermionOrderingSign(-1)*InternalFermionLoopSign(-1) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ | e- a d g QED^2==4 [{{2}} QCD] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"2 | (AutG(1))^(-1)*AntiFermionSpinSumSign(-1)*ExternalFermionOrderingSign(-1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(0,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(-1)*ExternalFermionOrderingSign(-1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(1,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(-1)*ExternalFermionOrderingSign(-1)*InternalFermionLoopSign(-1)) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e+ e- > d d~ / z QED==2 [{1}] --filter-selfenergies true",false)?,@"1 | (AutG(1))^(-1)*ExternalFermionOrderingSign(1) = 1");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "e+ e- > d d~ / z QED==2 [{1}] --filter-selfenergies false",false)?,@"3 | (AutG(1))^(-1)*3*ExternalFermionOrderingSign(1) = 3");

    // assert_snapshot!(feyngen_str(&mut cli, "amp", "e+ e- > d d~ / z | e- a d g QED^2==4 [{{1}} QCD] --compare-canonized-numerator true --number-of-samples-for-numerator-comparisons 3",false)?,@"10 | 4*UFO::GC_11^2*UFO::GC_1^(-2)+6");

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

    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --filter-zero-flow-edges false --fully-numerical-substitution-when-comparing-numerators false --compare-canonized-numerator true",false)?,@"10 | (AutG(1))^(-1)*7*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(10,24*G^2*TR*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(11,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(12,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(5,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(6,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(7,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(8,24*G^2*TR*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(9,12*G^2*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -12+-12*G^2*ee^(-2)+-48*G^2*TR*ee^(-2)");
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
    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > {d d~, g g} [{{1}}] | d g a --numerator-grouping only_detect_zeroes",false)?,@"1 | (AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1) = -1");
    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ g [{{2}}] | d g a --numerator-grouping only_detect_zeroes",false)?,@"3 | (AutG(1))^(-1)*3*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1) = -3");
    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ g [{{2}}] | d g a --symmetrize-left-right-states true",false)?,@"2 | (AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(0,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(1,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -3");
    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ z [{{2}}] | d g a --numerator-grouping only_detect_zeroes",false)?,@"0 | 0 = 0");


    // Full particle contents
    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{1}}] --symmetrize-left-right-states true",false)?,@"1 | (AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1) = -1");
    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states false --compare-canonized-numerator true --number-of-samples-for-numerator-comparisons 0 --fully-numerical-substitution-when-comparing-numerators false",true)?,@"10 | -12+-27*UFO::G^2*spenso::Nc*spenso::TR*UFO::ee^(-2)+27*UFO::G^2*spenso::Nc^-1*spenso::TR*UFO::ee^(-2)");
    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --compare-canonized-numerator true --number-of-samples-for-numerator-comparisons 0 --fully-numerical-substitution-when-comparing-numerators true",false)?,@"10 | -12+-36*UFO::G^2*UFO::ee^(-2)");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --compare-canonized-numerator false --number-of-samples-for-numerator-comparisons 3 --fully-numerical-substitution-when-comparing-numerators true",false)?,@"10 | -12+-36*UFO::G^2*UFO::ee^(-2)");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --no-compare-canonized-numerator --number-of-samples-for-numerator-comparisons 3 --fully-numerical-substitution-when-comparing-numerators false",false)?,@"10 | -12+-3*spenso::Nc*spenso::TR+3*spenso::Nc^-1*spenso::TR");

    // Only 1-flavour pure QCD corrections
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{1}}] --symmetrize-left-right-states true",false)?,@"1 | -1");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{2}} QCD=1] --symmetrize-left-right-states true",false)?,@"2 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{3}} QCD=2] --symmetrize-left-right-states true",false)?,@"16 | -41/2");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{4}} QCD=3] --symmetrize-left-right-states true",false)?,@"166 | -3107/14");


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
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{3}}] --numerator-grouping only_detect_zeroes",false)?,@"1321 | -1103/2");

    // Only 1-flavour pure QCD corrections
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{5}}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"6303 | -51683/24");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_sm_full_a_ddx() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm-full.json"))?;

    // Full particle contents
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{1}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"1 | -1");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"47 | -47");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"45 | -47");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"40 | -47");

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_slow_generate_sm_full_a_ddx() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm-full.json"))?;

    // Full particle contents
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ QED^2==4 [{{3}}] --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"339 | -368");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{3}}] --numerator-grouping no_grouping",false)?,@"8549 | -9111/2");

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_generate_sm_h_n_j() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Full particle contents
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g | h g b t ghg [{1}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"2 | -4");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g | h g b t ghg [{2}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"30 | -74");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g | h g b t ghg [{1}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"2 | -4");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g | h g b t ghg [{2}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"40 | -74");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g g | h g b t ghg [{1}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"6 | -24");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g g | h g b t ghg [{2}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"188 | -618");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g g | h g b t ghg [{1}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"4 | -24");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g g | h g b t ghg [{2}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"132 | -618");

    // Cross-section at 3-loops
    assert_snapshot!(feyngen_str(&mut cli, "xs", "h > g g | h g b t ghg [{{3}}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"8 | 8");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "h > g g | h g b t ghg [{{3}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"3 | 8");

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_slow_generate_sm_h_n_j_cross_section() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "xs", "h > g g g | h g b t ghg [{{4}}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"72 | 88");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "h > g g g | h g b t ghg [{{4}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_sign",false)?,@"26 | 88");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_slow_generate_sm_h_n_j() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g | h g b t ghg [{3}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"1194 | -1260");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "h > g g g | h g b t ghg [{3}] --symmetrize-left-right-states true --numerator-grouping only_detect_zeroes",false)?,@"6946 | -12429");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_epem_a_qqh() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Cross sections
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"3 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"3 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"30 | -30");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] --numerator-grouping only_detect_zeroes",false)?,@"594 | -261");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"21 | -21");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=2] --numerator-grouping only_detect_zeroes",false)?,@"171 | -171");

    // Enable grouping
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"16 | -30");
    // This grouping involves 21 new cancellations between graphs, so that modifies the total overall factor
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"266 | -327");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"20 | -21");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=1] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"157 | -171");

    // Enable left-right-symmetry
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > d d~ | d a e- ghg g QED^2==4 [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{3}} QCD=1] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"11 | -30");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g QED^2==6 [{{4}} QCD=2] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"166 | -327");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"16 | -21");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "e+ e- > b b~ h | d b h a e- ghg g z QED^2==6 [{{3}} QCD=1] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"104 | -171");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_a_qqh() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Cross sections
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d a ghg g QED^2==2 [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"3 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{2}}] --numerator-grouping only_detect_zeroes",false)?,@"3 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{3}} QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"30 | -30");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{4}} QCD=2] --numerator-grouping only_detect_zeroes",false)?,@"594 | -261");

    // Enable left-right-symmetry and grouping
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d a ghg g QED^2==2 [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{2}}] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{3}} QCD=1] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"11 | -30");
    // Additional cancellations modify the total overall factor
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | d b h a ghg g QED^2==4 [{{4}} QCD=2] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"166 | -327");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > b b~ h | b h a ghg g QED^2==4 [{{4}} QCD=2] --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"151 | -366");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_aa_ttx_amplitude() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{0} QCD=0] --numerator-grouping only_detect_zeroes",false)?,@"2 | 2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{1} QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"8 | 8");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{2} QCD=2] --numerator-grouping only_detect_zeroes",false)?,@"144 | 36");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{3} QCD=3] --numerator-grouping only_detect_zeroes",false)?,@"3424 | 28/3");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{0} QCD=0] --symmetrize-initial-states --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"1 | 2");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{1} QCD=1] --symmetrize-initial-states --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"4 | 8");
    // Six additional cancellations imply that the total overall factor is different
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{2} QCD=2] --symmetrize-initial-states --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"58 | 60");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "a a > t t~ | a t g b ghg QED==2 [{3} QCD=3] --symmetrize-initial-states --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"1221 | 1600/3");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_aa_ttx_cross_section() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --numerator-grouping only_detect_zeroes",false)?,@"4 | -4");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --numerator-grouping only_detect_zeroes",false)?,@"40 | -40");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] --numerator-grouping only_detect_zeroes",false)?,@"874 | -266");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --symmetrize-initial-states --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -4");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --symmetrize-initial-states --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"14 | -40");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] --symmetrize-initial-states --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"223 | -426");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --symmetrize-initial-states --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"2 | -4");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --symmetrize-initial-states --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"10 | -40");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g ghg QED^2==4 [{{3}} QCD=2] --symmetrize-initial-states --symmetrize-left-right-states true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling",false)?,@"128 | -502");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_generate_aa_ttx_cross_section_slow_filter() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{1}} QCD=0] --numerator-grouping only_detect_zeroes --max-multiplicity-for-fast-cut-filter 0",false)?,@"4 | -4");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{2}} QCD=1] --numerator-grouping only_detect_zeroes --max-multiplicity-for-fast-cut-filter 0",false)?,@"40 | -40");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a a > t t~ | a t g b ghg QED^2==4 [{{3}} QCD=2] --numerator-grouping only_detect_zeroes --max-multiplicity-for-fast-cut-filter 0",false)?,@"874 | -266");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_vacuum_amplitude_generation_kaapo() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm-default.json"))?;

    // 4-loop vaccuum contribution to the neutron start equation of state
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g d d~ ghG ghG~ [{4}] --numerator-grouping only_detect_zeroes --number-of-factorized-loop-subtopologies 1 1000 --number-of-fermion-loops 1 1000 --filter-snails false --filter-selfenergies false --filter-tadpoles false --max-n-bridges 0",false)?,@"53 | (AutG(1))^(-1)*2*ExternalFermionOrderingSign(1)+(AutG(1))^(-1)*7*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(12))^(-1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(2))^(-1)*15*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(2))^(-1)*2*ExternalFermionOrderingSign(1)+(AutG(3))^(-1)*4*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(4))^(-1)*3*ExternalFermionOrderingSign(1)+(AutG(4))^(-1)*9*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(6))^(-1)*2*ExternalFermionOrderingSign(1)+(AutG(6))^(-1)*2*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(8))^(-1)*5*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(8))^(-1)*ExternalFermionOrderingSign(1) = -179/12");

    Ok(())
}

#[test]
#[rustfmt::skip]
fn test_vacuum_amplitude_generation() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // Including all graphs
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{2}] --numerator-grouping only_detect_zeroes --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1",false)?,@"6 | -19/24");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1",false)?,@"5 | -19/24");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{3}] --numerator-grouping only_detect_zeroes --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1",false)?,@"36 | -83/48");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{3}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1",false)?,@"27 | -83/48");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{4}] --numerator-grouping only_detect_zeroes --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1",false)?,@"264 | -143/24");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{4}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges -1 --number-of-factorized-loop-subtopologies -1",false)?,@"179 | -731/168");

    // Only non-factorizable ones
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{2}] --numerator-grouping only_detect_zeroes --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"5 | -11/12");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"4 | -11/12");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{3}] --numerator-grouping only_detect_zeroes --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"29 | -23/16");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{3}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"21 | -23/16");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{4}] --numerator-grouping only_detect_zeroes --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"212 | -135/32");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{4}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"138 | -585/224");

    Ok(())
}

#[test]
#[rustfmt::skip]
#[ignore = "slow"]
fn test_slow_vacuum_amplitude_generation() -> Result<()> {
    let mut cli = get_test_cli(None, get_tests_workspace_path().join("feyn_gen_generation_test"), Some("feyngen".to_string()),true)?;
    cli.run_command(&format!("import model sm.json"))?;

    // CURRENTLY BUGGED BECAUSE OF SYMBOLICA INCORRECT SYMMETRY FACTORS. POSSIBLY OTHER SCENARIOS BUGGED TOO.
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{5}] --numerator-grouping only_detect_zeroes --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"2560 | -5785/384");
    assert_snapshot!(feyngen_str(&mut cli, "amp", "{} > {} | g ghg t u d QED==0 [{5}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1",false)?,@"1440 | 233015/51072");

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
