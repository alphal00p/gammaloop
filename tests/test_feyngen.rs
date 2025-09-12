use color_eyre::Result;

mod test_utils;
use gammalooprs::processes::ProcessCollection;
use symbolica::atom::Atom;
use symbolica::parse_lit;
use test_utils::{clean_test, get_test_cli, get_tests_workspace_path};

fn generate_graphs_and_count(process: String, model: String) -> Result<(usize, Atom)> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join("feyn_gen_generation_test"),
        Some("feyngen".to_string()),
    )?;
    cli.run_command(&format!("import model {}", model))?;
    cli.run_command(&format!(
        "generate graphs \"{}\" --clear-processes",
        process
    ))?;
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
    Ok((n_graphs, overall_factor_sum))
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
fn example_graph_count() -> Result<()> {
    assert_eq!(
        generate_graphs_and_count(
            "e+ e- > d d~ / Z QED^2==4 [{{1}} QCD]".to_string(),
            "sm.json".to_string(),
        )?,
        (2, parse_lit!(1))
    );
    Ok(())
}
