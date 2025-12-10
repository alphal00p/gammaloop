use bincode::de;
use color_eyre::Result;

mod test_utils;
use gammaloop_api::{commands::save, state::SyncSettings};
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

    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --filter-zero-flow-edges false --fully-numerical-substitution-when-comparing-numerators false --compare-canonized-numerator true",false)?,@"10 | (AutG(1))^(-1)*7*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(10,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(11,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(12,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(5,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(6,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(7,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(8,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(9,12*G^2*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -12+-12*G^2*ee^(-2)+-18*G^2*Nc*TR*ee^(-2)+18*G^2*Nc^(-1)*TR*ee^(-2)");
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

    // assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --numerator-grouping no_grouping --select-graphs GL03 GL04",true)?,@"2 | (AutG(1))^(-1)*2*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1) = -2");
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
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > {d d~, g g} [{{1}}] | d g a --numerator-grouping only_detect_zeroes",false)?,@"1 | (AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1) = -1");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ g [{{2}}] | d g a --numerator-grouping only_detect_zeroes",false)?,@"3 | (AutG(1))^(-1)*3*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ g [{{2}}] | d g a --symmetrize-left-right-states true",false)?,@"2 | (AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(0,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(1,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ z [{{2}}] | d g a --numerator-grouping only_detect_zeroes",false)?,@"0 | 0 = 0");


    // Full particle contents
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{1}}] --symmetrize-left-right-states true",false)?,@"1 | (AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1) = -1");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states false --compare-canonized-numerator true --number-of-samples-for-numerator-comparisons 0 --fully-numerical-substitution-when-comparing-numerators false",false)?,@"12 | (AutG(1))^(-1)*9*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(10,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(5,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(6,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(7,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(8,-9*G^2*Nc^(-1)*TR*ee^(-2)+9*G^2*Nc*TR*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(9,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -12+-27*G^2*Nc*TR*ee^(-2)+27*G^2*Nc^(-1)*TR*ee^(-2)");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --compare-canonized-numerator true --number-of-samples-for-numerator-comparisons 0 --fully-numerical-substitution-when-comparing-numerators true",false)?,@"12 | (AutG(1))^(-1)*9*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(10,12*G^2*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(5,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(6,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(7,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(8,12*G^2*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(9,12*G^2*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -12+-36*G^2*ee^(-2)");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --compare-canonized-numerator false --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --number-of-samples-for-numerator-comparisons 3 --fully-numerical-substitution-when-comparing-numerators false --symmetric-left-right-polarizations true",false)?,@"10 | (AutG(1))^(-1)*7*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(10,12*G^2*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(11,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(12,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(5,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(6,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(7,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(8,12*G^2*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(9,12*G^2*ee^(-2),(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -12+-36*G^2*ee^(-2)");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ [{{2}}] --symmetrize-left-right-states true --compare-canonized-numerator false --number-of-samples-for-numerator-comparisons 3 --fully-numerical-substitution-when-comparing-numerators true",false)?,@"13 | (AutG(1))^(-1)*11*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(5,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(6,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(8,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(9,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -15");

    // Only 1-flavour pure QCD corrections
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{1}}] --symmetrize-left-right-states true --symmetric-left-right-polarizations true",false)?,@"1 | (AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1) = -1");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{2}} QCD=1] --symmetrize-left-right-states true --symmetric-left-right-polarizations true",false)?,@"2 | (AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(0,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(1,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -3");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{3}} QCD=2] --symmetrize-left-right-states true --symmetric-left-right-polarizations true",false)?,@"16 | (AutG(1))^(-1)*4*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)+(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(0,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(1,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(10,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(11,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(12,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(13,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(14,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(15,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(16,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(17,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(18,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(19,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(23,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(24,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(3,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(4,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(5,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(6,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(8,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(9,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -41/2");
    assert_snapshot!(feyngen_str(&mut cli, "xs", "a > d d~ | d g ghG a QED^2==2 [{{4}} QCD=3] --symmetrize-left-right-states true --symmetric-left-right-polarizations true",true)?,@"158 | (AutG(1))^(-1)*10*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(1))^(-1)*3*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)+(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2)+(AutG(2))^(-1)*2*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(4))^(-1)*2*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+(AutG(6))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)+NumeratorDependentGrouping(0,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(1,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(10,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(104,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(105,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(106,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(107,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(108,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(11,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(111,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(112,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(115,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(116,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(118,1,(AutG(4))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(119,1,(AutG(4))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(12,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(120,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(121,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(122,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(123,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(124,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(125,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(126,1,(AutG(4))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(127,1,(AutG(4))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(128,1,(AutG(6))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(129,1,(AutG(6))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(13,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(130,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(131,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(132,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(133,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(134,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(135,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(136,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(137,2/7,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(138,2/7,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(139,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(14,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(140,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(141,2/7,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(142,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(143,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(144,2/7,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(145,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(146,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(147,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(148,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(149,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(15,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(150,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(151,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(152,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(153,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(154,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(155,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(156,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(157,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(158,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(159,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(16,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(160,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(161,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(162,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(163,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(164,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(165,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(166,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(167,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(168,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(169,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(17,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(172,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(173,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(174,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(175,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(176,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(177,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(178,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(179,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(18,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(180,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(181,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(182,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(183,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(184,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(185,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(186,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(187,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(188,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(189,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(19,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(190,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(191,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(192,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(193,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(194,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(195,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(196,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(197,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(198,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(199,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(2,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(20,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(200,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(201,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(202,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(203,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(204,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(205,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(206,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(207,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(208,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(209,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(21,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(215,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(216,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(217,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(218,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(22,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(221,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(222,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(223,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(224,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(227,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(228,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(229,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(23,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(230,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(231,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(232,-7/2,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(233,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(234,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(235,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(236,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(237,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(238,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(239,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(24,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(240,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(241,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(242,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(243,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(244,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(245,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(246,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(25,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(251,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(252,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(253,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(254,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(259,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(26,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(260,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(261,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(262,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(263,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(264,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(265,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(266,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(267,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(268,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(27,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(273,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(274,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(279,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(28,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(280,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(282,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(283,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(284,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(285,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(286,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(287,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(288,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(289,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(29,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(290,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(291,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(292,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(293,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(294,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(295,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(296,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(297,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(298,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(299,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(3,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(30,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(300,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(301,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(306,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(307,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(308,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(309,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(31,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(311,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(312,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(313,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(314,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(315,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(316,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(317,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(318,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(319,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(32,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(320,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(321,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(322,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(323,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(324,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(325,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(326,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(33,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(331,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(332,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(333,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(334,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(335,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(336,-2/7,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(337,-2/7,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(338,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(339,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(340,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(341,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(342,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(343,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(344,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(350,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(351,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(352,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(353,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(354,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(355,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(360,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(361,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(362,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(363,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(38,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(39,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(4,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(40,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(41,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(42,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(43,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(5,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(57,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(58,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(59,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(6,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(60,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(61,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(62,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(63,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(64,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(65,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(66,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(68,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(69,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(7,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(71,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(72,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(73,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(74,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(75,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(76,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(77,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(78,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(79,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(8,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(80,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(81,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(82,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(83,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(84,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(85,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(86,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(87,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(88,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(89,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(9,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*NumeratorIndependentSymmetryGrouping(2))+NumeratorDependentGrouping(90,1,(AutG(2))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(91,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(92,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(93,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(94,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(95,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(96,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1))+NumeratorDependentGrouping(97,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(98,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1))+NumeratorDependentGrouping(99,1,(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)) = -2939/14");


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
