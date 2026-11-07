use super::*;
use color_eyre::eyre::{ensure, eyre};
use gammaloop_api::{StateLoadOption, state::RunHistory};
use gammaloop_integration_tests::workspace_root;
use gammalooprs::graph::Graph;
use gammalooprs::processes::{CrossSection, ProcessCollection};
use gammalooprs::settings::global::GenerationSettings;
use gammalooprs::{settings::RuntimeSettings, utils::serde_utils::SmartSerde};
use linnet::half_edge::involution::HedgePair;
use std::fs;
use std::io::Result as IoResult;

fn is_generated_state_dir(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.contains("state"))
}

fn is_ignored_old_card(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.starts_with("OLD_"))
}

fn discover_toml_cards(root: &Path, cards: &mut Vec<PathBuf>) -> Result<()> {
    for entry in fs::read_dir(root)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            if !is_generated_state_dir(&path) {
                discover_toml_cards(&path, cards)?;
            }
        } else if !is_ignored_old_card(&path)
            // Drawing examples also contain Typst package manifests.
            && path.file_name().and_then(|name| name.to_str()) != Some("typst.toml")
            && path.extension().and_then(|extension| extension.to_str()) == Some("toml")
        {
            cards.push(path);
        }
    }
    Ok(())
}

#[test]
fn all_example_toml_cards_are_loadable() -> Result<()> {
    let workspace_root = workspace_root();
    let mut cards = Vec::new();
    for examples_dir in ["examples/api", "examples/cli"] {
        discover_toml_cards(&workspace_root.join(examples_dir), &mut cards)?;
    }
    cards.sort();
    if cards.is_empty() {
        return Err(eyre::eyre!("No example TOML cards were discovered"));
    }

    let failures = cards
        .iter()
        .filter_map(|card| {
            RunHistory::load(card)
                .map(|_| ())
                .or_else(|run_error| {
                    // Both schemas reject unknown fields, so a malformed run
                    // card cannot pass as a standalone runtime overlay.
                    RuntimeSettings::from_file(card, "example runtime settings")
                        .map(|_| ())
                        .map_err(|runtime_error| {
                            eyre::eyre!(
                                "Neither a valid run card nor runtime settings:\nRun card: {run_error:?}\nRuntime settings: {runtime_error:?}"
                            )
                        })
                })
                .err()
                .map(|err| {
                    let display_path = card.strip_prefix(&workspace_root).unwrap_or(card);
                    format!("{}:\n{err:?}", display_path.display())
                })
        })
        .collect_vec();

    if !failures.is_empty() {
        return Err(eyre::eyre!(
            "Failed to load {} example TOML card(s):\n\n{}",
            failures.len(),
            failures.join("\n\n")
        ));
    }
    Ok(())
}

/// The shipped tth@NNLO cross-section graphs import with the process
/// specification of their example cards, and their declared final states keep
/// every cut that this process accepts.
#[test]
fn epem_a_ttxh_nnlo_example_graphs_keep_every_process_cut() -> Result<()> {
    const SPEC: &str =
        "e+ e- > t t~ h | e+ e- g t t~ h d d~ ghG ghG~ a QCD^2==4 QED^2==6 [{{4}} QCD=2]";
    let nnlo = workspace_root().join("examples/cli/epem_a_ttxh/NNLO");
    assert!(
        fs::read_to_string(nnlo.join("gl297_optimized_lmbs.toml"))?
            .contains(&format!("--process-spec '{SPEC}'")),
        "the example cards must import the NNLO graphs with this process specification"
    );
    let state_path = get_tests_workspace_path().join("epem_a_ttxh_nnlo_example_graphs");
    let mut cli = get_test_cli(None, &state_path, None, true)?;
    cli.run_command("import model sm-default.json")?;
    cli.run_command(&format!(
        "import graphs {} --process-spec '{SPEC}' -p epem_a_tth -i NNLO",
        nnlo.join("graphs/GL297.dot").display()
    ))?;
    let model = cli.state.model.clone();
    let definition = cli.state.process_list.processes[0].definition.clone();

    // The Born final state with up to two further particles of the
    // specification: a superset of every final state the process accepts.
    let extra = ["g", "d", "d~", "ghG", "ghG~", "a"];
    let mut every_final_state = vec!["H,t,t~".to_string()];
    for (position, first) in extra.iter().enumerate() {
        every_final_state.push(format!("H,t,t~,{first}"));
        every_final_state.extend(
            extra[position..]
                .iter()
                .map(|second| format!("H,t,t~,{first},{second}")),
        );
    }
    let every_final_state = format!("final_state = \"{}\";", every_final_state.join("|"));
    // (cuts finalized at import, cuts compatible with the process)
    let cut_counts = |source: &str| -> Result<(usize, usize)> {
        let graphs = Graph::import_dot(source, &model)?;
        let imported = graphs[0].finalized_cuts.len();
        let cross_section = CrossSection::from_graph_list("check".into(), graphs)?;
        let count = cross_section.supergraphs[0].cutkosky_cut_count_for_process(
            &model,
            &definition,
            &GenerationSettings::default(),
        )?;
        Ok((imported, count.process_compatible_cuts))
    };
    let check = |path: &Path| -> Result<()> {
        let source = fs::read_to_string(path)?;
        let declared = source
            .lines()
            .map(str::trim)
            .find(|line| line.starts_with("final_state"))
            .ok_or_else(|| eyre!("no final_state attribute"))?;
        let shipped = cut_counts(&source)?;
        let reference = cut_counts(&source.replacen(declared, &every_final_state, 1))?;
        ensure!(
            shipped == reference,
            "declared final states keep {shipped:?} (imported, process-compatible) cuts, \
             every final state keeps {reference:?}"
        );
        Ok(())
    };

    let mut paths = fs::read_dir(nnlo.join("graphs"))?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<IoResult<Vec<_>>>()?;
    paths.sort();
    assert_eq!(paths.len(), 166);
    let failures = paths
        .iter()
        .filter_map(|path| {
            check(path)
                .err()
                .map(|error| format!("{}: {error:?}", path.display()))
        })
        .collect_vec();
    clean_test(&state_path);

    ensure!(
        failures.is_empty(),
        "{} NNLO example graph(s) failed:\n\n{}",
        failures.len(),
        failures.join("\n\n")
    );
    Ok(())
}

#[test]
fn test_scalar_bubble_example_cli() -> Result<()> {
    let state_path = get_tests_workspace_path().join("scalar_bubble_example");
    let mut cli = get_example_cli(
        "scalar_topologies/bubble.toml",
        &["generate"],
        Some(state_path.clone()),
        None,
        true,
    )?;

    assert!(
        !cli.state.process_list.processes.is_empty(),
        "No processes were generated"
    );
    assert_eq!(
        cli.state.model.name(),
        "scalars",
        "Expected scalars model to be loaded"
    );
    assert!(
        state_path.join("state_manifest.toml").exists(),
        "Expected saved state manifest in {}",
        state_path.display()
    );

    let inspections = ["bubble", "bubble_no_integrated_UV"].map(|process| Inspect {
        process: Some(ProcessRef::Unqualified(process.to_string())),
        integrand_name: Some("scalar_bubble_below_thres".to_string()),
        point: vec![0.11, 0.23, 0.37],
        discrete_dim: vec![0, 0],
        ..Default::default()
    });
    let expected = inspections
        .iter()
        .map(|inspect| inspect.run(&mut cli.state))
        .collect::<Result<Vec<_>>>()?;
    let mut reloaded = StateLoadOption::read_only(&state_path).load()?;
    for (inspect, expected) in inspections.iter().zip(expected) {
        assert!(
            expected.1.re.is_finite()
                && expected.1.im.is_finite()
                && expected.1.re.hypot(expected.1.im) > 0.0,
            "state roundtrip requires a finite, nonzero reference"
        );
        assert_eq!(
            inspect.run(&mut reloaded.state)?,
            expected,
            "saved amplitude changed after reload"
        );
    }

    clean_test(&state_path);
    Ok(())
}

/// The shipped subtracted qq -> aaa graphs import with the model and options of
/// their example card, with all 19 graphs that GammaLoop imported before FeynKit.
#[test]
fn qq_aaa_subtracted_example_card_imports_all_graphs() -> Result<()> {
    const IMPORT: &str =
        "import graphs ./examples/cli/qq_aaa/qqx_aaa_subtracted.dot -p qq_aaa -i subtracted";
    let card = example_run_card("qq_aaa/qqx_aaa_subtracted.toml")?;
    let load = card
        .command_blocks
        .iter()
        .find(|block| block.name == "load")
        .expect("the qq_aaa card has a load block");
    let commands = load
        .commands
        .iter()
        .filter_map(|command| command.raw_string.as_deref())
        .collect_vec();
    assert!(
        commands.contains(&"import model sm-default.json") && commands.contains(&IMPORT),
        "the example card must load these graphs under this model and these options"
    );

    let state_path = get_tests_workspace_path().join("qq_aaa_subtracted_example_import");
    let mut cli = get_test_cli(None, &state_path, None, true)?;
    run_commands(&mut cli, &["import model sm-default.json", IMPORT])?;
    let ProcessCollection::Amplitudes(amplitudes) = &cli.state.process_list.processes[0].collection
    else {
        panic!("the qq_aaa example imports an amplitude")
    };
    assert_eq!(amplitudes["subtracted"].graphs.len(), 19);

    clean_test(&state_path);
    Ok(())
}

/// The one-loop photon self-energy card divides its transverse projector by
/// `Q(0,mink(4,1))*Q(0,mink(4,1))`, which admission writes as the even power
/// `Q(0,mink(4,1))^-2`. Its value matches the same projector with `Q^2`
/// written by components, and the values from before the typed rewrite.
#[test]
fn a_a_one_loop_projector_with_contracted_inverse_square_matches_components() -> Result<()> {
    const INVARIANT: &str =
        "gammalooprs::Q(0,spenso::mink(4,1))*gammalooprs::Q(0,spenso::mink(4,1))";
    const COMPONENTS: &str = "(gammalooprs::Q(0,spenso::cind(0))^2-gammalooprs::Q(0,spenso::cind(1))^2\
        -gammalooprs::Q(0,spenso::cind(2))^2-gammalooprs::Q(0,spenso::cind(3))^2)";
    let card = example_run_card("a_a/1L/a_a.toml")?;
    let block = |name: &str| {
        card.command_blocks
            .iter()
            .find(|block| block.name == name)
            .unwrap_or_else(|| panic!("the a_a card has a {name} block"))
            .commands
            .iter()
            .filter_map(|command| command.raw_string.clone())
            .filter(|command| !command.starts_with("save "))
            .collect_vec()
    };
    let generate = block("generate_diagrams");
    let shipped = generate
        .iter()
        .find(|command| command.contains(INVARIANT) && command.contains("-p ff_pi -i 1L"))
        .expect("the a_a card generates ff_pi with the contracted inverse square");
    let components = shipped
        .replace(INVARIANT, COMPONENTS)
        .replace("-p ff_pi -i 1L", "-p ff_pi_components -i 1L");
    assert_eq!(shipped.matches(INVARIANT).count(), 2);

    let state_path = get_tests_workspace_path().join("a_a_one_loop_inverse_square_projector");
    // The card's own generation and runtime settings, without running its blocks.
    let mut cli = get_example_cli("a_a/1L/a_a.toml", &[], Some(state_path.clone()), None, true)?;
    let commands = generate
        .iter()
        .cloned()
        .chain([components])
        .chain(block("generate_integrands"))
        .chain(["generate existing -p ff_pi_components -i 1L".to_owned()])
        .chain(block("set_model_parameters"))
        .chain(block("set_kinematics"))
        .collect_vec();
    run_commands(&mut cli, &commands.iter().map(String::as_str).collect_vec())?;

    let inspect = |process: &str, momentum_space: bool| Inspect {
        process: Some(ProcessRef::Unqualified(process.to_owned())),
        integrand_name: Some("1L".to_owned()),
        point: vec![0.31, 0.17, 0.43],
        momentum_space,
        graph_id: momentum_space.then_some(0),
        discrete_dim: if momentum_space { vec![] } else { vec![0] },
        ..Default::default()
    };
    // Values of `inspect -p name:ff_pi -i 1L` with the main-branch CLI.
    for (momentum_space, before_rewrite) in
        [(true, 1.049327628450592e-3), (false, 8.351826969293263e-3)]
    {
        let (_, shipped) = inspect("ff_pi", momentum_space).run(&mut cli.state)?;
        let (_, components) = inspect("ff_pi_components", momentum_space).run(&mut cli.state)?;
        for value in [shipped, components] {
            assert!(
                (value.re - before_rewrite).abs() <= 1e-12 * before_rewrite && value.im == 0.0,
                "momentum_space={momentum_space}: {value} differs from {before_rewrite}"
            );
        }
        assert!((shipped.re - components.re).abs() <= 1e-13 * before_rewrite);
    }

    clean_test(&state_path);
    Ok(())
}

/// Relation test for `Q3(e, mink(4, mu))`, the spatial part (0, q) of an edge momentum, in a
/// supplied projector. The one-loop `ff_pi_0` card contracts r = q - 2*q3 = (q0, -q), with q the
/// incoming photon momentum. In one build and at one point:
/// - r written without `Q3`, as 2*q0*delta(cind(0), mu) - q, gives the same value;
/// - r with the spatial part negated, 3*q - 2*q0*delta(cind(0), mu) = (q0, 3q), does not;
/// - the incoming and outgoing photon energies `Q(0, cind(0))` and `Q(1, cind(0))` are equal.
#[test]
fn a_a_one_loop_spatial_momentum_projector_matches_components() -> Result<()> {
    const PROJECTOR_FLAG: &str = "--global-prefactor-projector '";
    let card = example_run_card("a_a/1L/a_a.toml")?;
    let block = |name: &str| {
        card.command_blocks
            .iter()
            .find(|block| block.name == name)
            .unwrap_or_else(|| panic!("the a_a card has a {name} block"))
            .commands
            .iter()
            .filter_map(|command| command.raw_string.clone())
            .filter(|command| !command.starts_with("save "))
            .collect_vec()
    };
    let generate = block("generate_diagrams");
    let shipped = generate
        .iter()
        .find(|command| command.contains("-p ff_pi_0 -i 1L"))
        .expect("the a_a card generates ff_pi_0");
    let projector_start = shipped
        .find(PROJECTOR_FLAG)
        .expect("ff_pi_0 has a supplied projector")
        + PROJECTOR_FLAG.len();
    let projector_end = projector_start + shipped[projector_start..].find('\'').unwrap();
    let projector = &shipped[projector_start..projector_end];
    assert!(projector.contains("gammalooprs::Q3("));
    // The photon Lorentz slots, in the spelling this card uses.
    let slots = projector
        .match_indices("spenso::mink(4,gammalooprs::hedge(")
        .map(|(start, _)| {
            let slot = &projector[start + "spenso::mink(4,".len()..];
            &slot[..slot.find(')').unwrap() + 1]
        })
        .unique()
        .collect_vec();
    assert_eq!(slots.len(), 2, "two photon slots in {projector}");

    let q0 = "gammalooprs::Q(0,spenso::cind(0))";
    let q = |slot: &str| format!("gammalooprs::Q(0,spenso::mink(4,{slot}))");
    let temporal = |slot: &str| format!("{q0}*spenso::δ(spenso::cind(0),spenso::mink(4,{slot}))");
    let vector = |name: &str, slot: &str| match name {
        "ff_pi_0_spatial" => format!(
            "({} - 2*gammalooprs::Q3(0,spenso::mink(4,{slot})))",
            q(slot)
        ),
        "ff_pi_0_temporal" => format!("(2*{} - {})", temporal(slot), q(slot)),
        "ff_pi_0_negated" => format!("(3*{} - 2*{})", q(slot), temporal(slot)),
        _ => unreachable!("unknown projector vector {name}"),
    };
    let mut processes = ["ff_pi_0_spatial", "ff_pi_0_temporal", "ff_pi_0_negated"]
        .map(|name| {
            let [a, b] = [slots[0], slots[1]].map(|slot| vector(name, slot));
            (name, format!("(-1/(4*{q0}^4))*{a}*{b}"))
        })
        .to_vec();
    let temporal_projector = processes[1].1.clone();
    for (name, energy) in [
        ("ff_pi_0_incoming_energy", q0.to_owned()),
        (
            "ff_pi_0_outgoing_energy",
            "gammalooprs::Q(1,spenso::cind(0))".to_owned(),
        ),
    ] {
        processes.push((name, format!("{energy}*{temporal_projector}")));
    }

    let state_path = get_tests_workspace_path().join("a_a_one_loop_spatial_momentum_projector");
    // The card's own generation and runtime settings, without running its blocks.
    let mut cli = get_example_cli("a_a/1L/a_a.toml", &[], Some(state_path.clone()), None, true)?;
    let commands = generate
        .iter()
        .filter(|command| !command.contains("generate amp"))
        .cloned()
        .chain(processes.iter().map(|(name, projector)| {
            format!(
                "{}{projector}{}",
                &shipped[..projector_start],
                &shipped[projector_end..]
            )
            .replace("-p ff_pi_0 -i 1L", &format!("-p {name} -i 1L"))
        }))
        .chain(
            processes
                .iter()
                .map(|(name, _)| format!("generate existing -p {name} -i 1L")),
        )
        .chain(block("set_model_parameters"))
        .chain(block("set_kinematics_onshell"))
        .collect_vec();
    run_commands(&mut cli, &commands.iter().map(String::as_str).collect_vec())?;

    for momentum_space in [true, false] {
        let mut value = |process: &str| -> Result<(f64, f64)> {
            let inspect = Inspect {
                process: Some(ProcessRef::Unqualified(process.to_owned())),
                integrand_name: Some("1L".to_owned()),
                point: vec![0.31, 0.17, 0.43],
                momentum_space,
                graph_id: momentum_space.then_some(0),
                discrete_dim: if momentum_space { vec![] } else { vec![0] },
                ..Default::default()
            };
            let (_, value) = inspect.run(&mut cli.state)?;
            Ok((value.re, value.im))
        };
        let distance = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0).hypot(a.1 - b.1);
        let spatial = value("ff_pi_0_spatial")?;
        let temporal = value("ff_pi_0_temporal")?;
        let negated = value("ff_pi_0_negated")?;
        let incoming = value("ff_pi_0_incoming_energy")?;
        let outgoing = value("ff_pi_0_outgoing_energy")?;
        let scale = temporal.0.hypot(temporal.1);
        assert!(scale > 0.0, "momentum_space={momentum_space}");
        assert!(
            distance(spatial, temporal) <= 1e-12 * scale,
            "momentum_space={momentum_space}: Q3 form {spatial:?} differs from {temporal:?}"
        );
        assert!(
            distance(negated, temporal) >= 1e-2 * scale,
            "momentum_space={momentum_space}: negated form {negated:?} matches {temporal:?}"
        );
        assert!(
            outgoing.0.hypot(outgoing.1) > 0.0
                && distance(incoming, outgoing) <= 1e-12 * outgoing.0.hypot(outgoing.1),
            "momentum_space={momentum_space}: incoming energy {incoming:?} differs from outgoing \
             {outgoing:?}"
        );
    }

    clean_test(&state_path);
    Ok(())
}

/// Relation test for an outgoing leg's `Q3` in a cross-section projector. In
/// e+ e- -> d d~ at LO the outgoing copies 2 and 3 of the initial states are sewn
/// to the incoming legs 0 and 1, so `Q(2, ...)` and `Q3(2, ...)` are the e+ momentum
/// and its spatial part. In one build and at one point, with the massless spin-sum
/// projector times a factor:
/// - `Q3(2, mu) * delta(cind(3), mu)` equals its explicit-component twin
///   `(Q(2, mu) - Q(2, cind(0)) * delta(cind(0), mu)) * delta(cind(3), mu)`;
/// - it equals the incoming leg's `Q3(0, mu) * delta(cind(3), mu)`, and
///   `Q3(2, mu) * Q(2, mu)` equals `Q3(0, mu) * Q(0, mu)`;
/// - the twin with the spatial part negated differs.
#[test]
fn epem_a_ddx_cross_section_outgoing_leg_spatial_momentum_projector_matches_components()
-> Result<()> {
    const PROCESS_FLAG: &str = "-p epem_a_ddx -i LO";
    let card = example_run_card("epem_a_ddx/LO/epem_a_ddx.toml")?;
    let block = |name: &str| {
        card.command_blocks
            .iter()
            .find(|block| block.name == name)
            .unwrap_or_else(|| panic!("the epem_a_ddx card has a {name} block"))
            .commands
            .iter()
            .filter_map(|command| command.raw_string.clone())
            .filter(|command| !command.starts_with("save "))
            .collect_vec()
    };
    let generate = block("generate");
    let shipped = generate
        .iter()
        .find(|command| command.starts_with("generate xs") && command.contains(PROCESS_FLAG))
        .expect("the epem_a_ddx card generates the LO cross section");

    // Process legs: e+ = 0, e- = 1, and their sewn outgoing copies 2 and 3. The
    // massless spin sums are u(p) ubar(p) and v(p) vbar(p) closed over each pair.
    let spin_sum = String::from(
        "1/4*gammalooprs::Q(0,spenso::mink(4,101))*gammalooprs::Q(1,spenso::mink(4,102))\
         *spenso::gamma(spenso::bis(4,gammalooprs::hedge(1)),spenso::bis(4,gammalooprs::hedge(3)),\
         spenso::mink(4,102))\
         *spenso::gamma(spenso::bis(4,gammalooprs::hedge(2)),spenso::bis(4,gammalooprs::hedge(0)),\
         spenso::mink(4,101))",
    );
    let mu = "spenso::mink(4,7)";
    let z = format!("spenso::δ(spenso::cind(3),{mu})");
    let q = |leg: usize| format!("gammalooprs::Q({leg},{mu})");
    let q0 = |leg: usize| format!("gammalooprs::Q({leg},spenso::cind(0))");
    let q3 = |leg: usize| format!("gammalooprs::Q3({leg},{mu})");
    let t = format!("spenso::δ(spenso::cind(0),{mu})");
    let factors = [
        ("xs_spatial_out_z", format!("{}*{z}", q3(2))),
        (
            "xs_temporal_out_z",
            format!("({} - {}*{t})*{z}", q(2), q0(2)),
        ),
        (
            "xs_negated_out_z",
            format!("({}*{t} - {})*{z}", q0(2), q(2)),
        ),
        ("xs_spatial_in_z", format!("{}*{z}", q3(0))),
        ("xs_spatial_out_square", format!("{}*{}", q3(2), q(2))),
        ("xs_spatial_in_square", format!("{}*{}", q3(0), q(0))),
    ];

    let state_path =
        get_tests_workspace_path().join("epem_a_ddx_cross_section_spatial_momentum_projector");
    // The card's own generation and runtime settings, without running its blocks.
    let mut cli = get_example_cli(
        "epem_a_ddx/LO/epem_a_ddx.toml",
        &[],
        Some(state_path.clone()),
        None,
        true,
    )?;
    // The card's setup (model import) first, then its model point, which decides
    // the threshold counterterms at generation, as in its own generate block.
    let commands = generate
        .iter()
        .filter(|command| {
            !command.starts_with("generate ")
                && !command.starts_with("duplicate ")
                && !command.starts_with("run ")
        })
        .cloned()
        .chain(block("set_model_parameters"))
        .chain(factors.iter().map(|(name, factor)| {
            shipped
                .replace(
                    "--numerator-grouping group_identical_graphs_up_to_scalar_rescaling",
                    "--numerator-grouping no_grouping",
                )
                .replace(
                    PROCESS_FLAG,
                    &format!(
                        "--global-prefactor-projector '({factor})*{spin_sum}' -p {name} -i LO"
                    ),
                )
        }))
        .collect_vec();
    run_commands(&mut cli, &commands.iter().map(String::as_str).collect_vec())?;

    for momentum_space in [true, false] {
        let mut value = |process: &str| -> Result<(f64, f64)> {
            let inspect = Inspect {
                process: Some(ProcessRef::Unqualified(process.to_owned())),
                integrand_name: Some("LO".to_owned()),
                point: if momentum_space {
                    vec![130.0, -70.0, 210.0]
                } else {
                    vec![0.31, 0.17, 0.43]
                },
                momentum_space,
                graph_id: momentum_space.then_some(0),
                // One graph: the LO cross section has no discrete graph dimension.
                ..Default::default()
            };
            let (_, value) = inspect.run(&mut cli.state)?;
            Ok((value.re, value.im))
        };
        let distance = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0).hypot(a.1 - b.1);
        let close = |a: (f64, f64), b: (f64, f64)| distance(a, b) <= 1e-12 * b.0.hypot(b.1);
        let spatial_out = value("xs_spatial_out_z")?;
        let temporal_out = value("xs_temporal_out_z")?;
        let negated_out = value("xs_negated_out_z")?;
        let spatial_in = value("xs_spatial_in_z")?;
        let square_out = value("xs_spatial_out_square")?;
        let square_in = value("xs_spatial_in_square")?;
        let scale = temporal_out.0.hypot(temporal_out.1);
        assert!(scale > 0.0, "momentum_space={momentum_space}");
        assert!(
            close(spatial_out, temporal_out),
            "momentum_space={momentum_space}: Q3 form {spatial_out:?} differs from {temporal_out:?}"
        );
        assert!(
            distance(negated_out, temporal_out) >= 1e-2 * scale,
            "momentum_space={momentum_space}: negated form {negated_out:?} matches {temporal_out:?}"
        );
        assert!(
            close(spatial_out, spatial_in),
            "momentum_space={momentum_space}: outgoing Q3 {spatial_out:?} differs from incoming \
             {spatial_in:?}"
        );
        assert!(
            square_in.0.hypot(square_in.1) > 0.0 && close(square_out, square_in),
            "momentum_space={momentum_space}: outgoing Q3.Q {square_out:?} differs from incoming \
             {square_in:?}"
        );
    }

    clean_test(&state_path);
    Ok(())
}

#[test]
fn test_aa_aa_example_card() -> Result<()> {
    let run_history = example_run_card("aa_aa/1L/aa_aa.toml")?;

    assert_eq!(
        run_history.command_blocks[0].commands[0]
            .raw_string
            .as_deref(),
        Some("import model sm-default.json")
    );
    assert_eq!(
        run_history.cli_settings.state.folder,
        PathBuf::from("./examples/cli/aa_aa/1L/gammaloop_state")
    );
    assert_eq!(run_history.default_runtime_settings.kinematics.e_cm, 300.0);
    Ok(())
}

#[test]
fn test_epem_a_tth_nlo_example_cli() -> Result<()> {
    let state_path = get_tests_workspace_path().join("epem_a_tth_nlo_example");
    let cli = get_example_cli(
        "epem_a_ttxh/NLO/epem_a_tth_NLO.toml",
        &["generate_diagrams"],
        Some(state_path.clone()),
        None,
        true,
    )?;

    assert_eq!(
        cli.state.model.name(),
        "sm",
        "Expected SM model to be loaded"
    );
    assert!(
        !cli.state.process_list.processes.is_empty(),
        "No processes were generated"
    );
    assert_eq!(
        cli.state.process_list.processes[0].definition.folder_name,
        "epem_a_tth"
    );
    assert!(
        state_path.join("justfile").exists(),
        "Expected dot export helper files to be created in {}",
        state_path.display()
    );

    clean_test(&state_path);
    Ok(())
}

/// The explicit counterterm example keeps its nine external legs, including the
/// dummy momentum carriers used by its supplied numerator expressions.
#[test]
fn qq_hhh_two_loop_example_keeps_runtime_numerators_and_dummy_legs() -> Result<()> {
    const COMMANDS: [&str; 7] = [
        "import model sm-default.json",
        "set model MT=173.0",
        "set model MH=125.0",
        "set model WT=0.0",
        "set model WH=0.0",
        "set model ymt=173.0",
        "import graphs ./examples/cli/qq_hhh_2L/qq_hhh_2L.dot -p qq_hhh_2L -i subtracted",
    ];
    let card = example_run_card("qq_hhh_2L/qq_hhh_2L.toml")?;
    let load = card
        .command_blocks
        .iter()
        .find(|block| block.name == "load")
        .expect("the qq_hhh card has a load block");
    let commands = load
        .commands
        .iter()
        .filter_map(|command| command.raw_string.as_deref())
        .collect_vec();
    assert!(COMMANDS.iter().all(|command| commands.contains(command)));

    let state_path = get_tests_workspace_path().join("qq_hhh_two_loop_example_import");
    let mut cli = get_test_cli(None, &state_path, None, true)?;
    run_commands(&mut cli, &COMMANDS)?;
    let ProcessCollection::Amplitudes(amplitudes) = &cli.state.process_list.processes[0].collection
    else {
        panic!("the qq_hhh example imports an amplitude")
    };
    let graphs = &amplitudes["subtracted"].graphs;
    let names = graphs
        .iter()
        .map(|graph| graph.graph.name.as_str())
        .sorted()
        .collect_vec();
    assert_eq!(
        names,
        [
            "GL05",
            "GL05_isr_p1_ct",
            "GL05_isr_p2_ct",
            "GL06",
            "GL06_isr_p1_ct",
            "GL06_isr_p2_ct",
            "GL12",
            "GL12_isr_p1_ct",
            "GL12_isr_p2_ct",
            "GL13",
            "GL13_isr_p1_ct",
            "GL13_isr_p2_ct",
        ]
    );
    for graph in graphs {
        assert_eq!(graph.graph.iter_edges().count(), 17, "{}", graph.graph.name);
        assert_eq!(
            graph
                .graph
                .iter_edges()
                .filter(|(pair, _, _)| matches!(pair, HedgePair::Unpaired { .. }))
                .count(),
            9,
            "{}",
            graph.graph.name,
        );
        assert!(
            graph
                .graph
                .underlying
                .iter_nodes()
                .all(|(_, _, vertex)| vertex.vertex_rule.is_none())
        );
    }
    let dummy_count = graphs
        .iter()
        .flat_map(|graph| graph.graph.iter_edges())
        .filter(|(_, _, edge)| edge.data.is_dummy)
        .count();
    assert_eq!(dummy_count, 32);

    clean_test(&state_path);
    Ok(())
}
