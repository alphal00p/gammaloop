use color_eyre::Result;

mod test_utils;
use gammalooprs::processes::ProcessCollection;
use symbolica::atom::{Atom, AtomCore};
use symbolica::parse_lit;
use test_utils::{clean_test, get_test_cli, get_tests_workspace_path};

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
    (n_graphs, overall_factor_sum)
}

fn generate_graphs_and_count(
    generation_type: &str,
    process: &str,
    model: &str,
) -> Result<(usize, Atom)> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join("feyn_gen_generation_test"),
        Some("feyngen".to_string()),
    )?;
    cli.run_command(&format!("import model {}", model))?;
    cli.run_command(&format!(
        "generate {} \"{}\" --clear-existing-processes --only-diagrams",
        generation_type, process
    ))?;
    Ok(count_graphs_in_processes(&cli))
}

#[test]
fn simple_epem_ddx_generation() -> Result<()> {
    let mut cli = get_test_cli(
        Some("sm_load.toml".into()),
        get_tests_workspace_path().join("epem_ddx_generation"),
        Some("epem_ddx_generation".to_string()),
    )?;

    cli.run_command("generate amp e+ e- > d d~")?;

    clean_test(&cli);

    Ok(())
}

#[test]
fn graph_count_from_amplitude_load() -> Result<()> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join("feyn_gen_generation_test"),
        Some("feyngen".to_string()),
    )?;
    cli.run_command(&format!("import model sm.json"))?;
    cli.run_command(&format!(
        "import amplitude ./tests/resources/graphs/qqx_aaa_subtracted.dot",
    ))?;
    let (n_graphs, overall_factor_sum) = count_graphs_in_processes(&cli);

    assert_snapshot!(format!("{n_graphs}"),@"19");
    assert_snapshot!(overall_factor_sum.to_canonical_string(),@"8-11𝑖");

    Ok(())
}

#[test]
fn example_graph_count() -> Result<()> {
    let (n_graphs, overall_factor_sum) =
        generate_graphs_and_count("xs", "e+ e- > d d~ / Z QED^2==4 [{{1}} QCD]", "sm.json")?;
    assert_snapshot!(format!("{n_graphs}"),@"0");
    assert_snapshot!(overall_factor_sum.to_canonical_string(),@"expr");
    Ok(())
}
