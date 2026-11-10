//! Validation of the ordinary-R pole identity on Vakint's exact 1L-4L masters.
//!
//! `K G = -Delta(G) - K sum_S Delta(S) G/S` is compared with Vakint's direct
//! evaluation of the master graph in the same conventions. Selection and
//! output are controlled by environment variables (the tests are ignored by
//! default):
//! - `RSTAR_MASSIVE_MASTERS`: comma-separated `<L>L/<P>` or `<L>L/*`;
//! - `RSTAR_MASSIVE_OUT`: JSON-lines output file (appended);
//! - `RSTAR_MASSIVE_ROUTES`: `i`, `ii` or `i,ii` for Delta(G) (default `i,ii`);
//! - `RSTAR_MASSIVE_EXTRA`: extra epsilon terms for the quotients (default 0).

use std::{io::Write, path::PathBuf, time::Instant};

use eyre::{Result, eyre};
use linnet::half_edge::involution::EdgeIndex;
use symbolica::{
    atom::{Atom, AtomCore},
    symbol,
};
use vakint::Vakint;

use crate::{
    initialisation::test_initialise,
    utils::{GS, load_generic_model},
};

use super::{
    RStarPlan,
    massive::{VacuumTopology, massive_vacuum_direct},
};

fn vakint_data_version() -> String {
    let manifest: serde_json::Value = serde_json::from_slice(include_bytes!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../vakint/data/manifest.json"
    )))
    .expect("vakint's vendored data manifest is JSON");
    manifest["version"].as_str().unwrap().to_owned()
}

fn data_root() -> PathBuf {
    std::env::var_os("VAKINT_TEST_EXTERNAL_DATA")
        .map(PathBuf::from)
        .expect("VAKINT_TEST_EXTERNAL_DATA must name a local data release")
}

pub(super) fn loaded_vakint(max_loops: usize) -> Result<Vakint> {
    let mut vakint = Vakint::new()?;
    vakint.load_official_external_data(
        max_loops,
        &vakint::OfficialExternalDataOptions {
            version: vakint_data_version(),
            local_artifacts_path: Some(data_root()),
            ..Default::default()
        },
    )?;
    Ok(vakint)
}

/// Parse `vA -> vB [id=N, ... power=P ...]` lines of a master's graph.dot.
pub(super) fn master_topology(loops: usize, master: &str) -> Result<VacuumTopology> {
    let path = data_root().join(format!("{loops}L/masters/{master}/graph.dot"));
    let text = std::fs::read_to_string(&path)?;
    let mut edges = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        let Some((left, rest)) = line.split_once(" -> ") else {
            continue;
        };
        let Some((right, attributes)) = rest.split_once(" [") else {
            continue;
        };
        let vertex = |name: &str| -> Result<usize> {
            name.trim()
                .strip_prefix('v')
                .ok_or_else(|| eyre!("unexpected vertex {name}"))?
                .parse::<usize>()
                .map_err(Into::into)
        };
        let power = attributes
            .split([',', ' ', ']'])
            .find_map(|field| field.strip_prefix("power="))
            .ok_or_else(|| eyre!("no power in {line}"))?
            .parse::<u32>()?;
        edges.push((vertex(left)?, vertex(right)?, power));
    }
    let topology = VacuumTopology {
        name: format!("m{loops}L_{master}"),
        edges,
    };
    eyre::ensure!(topology.loops() == loops, "{master}: wrong loop count");
    Ok(topology)
}

fn selection() -> Vec<(usize, String)> {
    let spec = std::env::var("RSTAR_MASSIVE_MASTERS").unwrap_or_else(|_| "2L/P2".into());
    let mut out = Vec::new();
    for item in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        let (loops, master) = item.split_once("L/").expect("<L>L/<P>");
        let loops: usize = loops.parse().unwrap();
        if master == "*" {
            let mut names = std::fs::read_dir(data_root().join(format!("{loops}L/masters")))
                .unwrap()
                .map(|entry| entry.unwrap().file_name().into_string().unwrap())
                .collect::<Vec<_>>();
            names.sort_by_key(|name| name[1..].parse::<usize>().unwrap());
            out.extend(names.into_iter().map(|name| (loops, name)));
        } else {
            out.push((loops, master.to_owned()));
        }
    }
    out
}

fn unit_scales(expression: &Atom, mass: symbolica::atom::Symbol) -> Atom {
    expression
        .replace(mass)
        .with(Atom::one())
        .replace(GS.mu_r_sq)
        .with(Atom::one())
        .expand()
}

#[test]
#[ignore = "needs VAKINT_TEST_EXTERNAL_DATA; selection via RSTAR_MASSIVE_MASTERS"]
fn massive_vacuum_masters_reproduce_their_poles() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    let settings = super::tests::settings();
    let routes = std::env::var("RSTAR_MASSIVE_ROUTES").unwrap_or_else(|_| "i,ii".into());
    let extra: i64 = std::env::var("RSTAR_MASSIVE_EXTRA")
        .ok()
        .map(|v| v.parse().unwrap())
        .unwrap_or(0);
    let out = std::env::var_os("RSTAR_MASSIVE_OUT").map(PathBuf::from);
    let selected = selection();
    let max_loops = selected.iter().map(|(l, _)| *l).max().unwrap_or(2).max(3);
    let vakint = loaded_vakint(max_loops)?;
    let mass_symbol = symbol!("rstar_vacuum_mass"; Scalar);
    let mass = Atom::var(mass_symbol);
    let mut failures = Vec::new();
    for (loops, master) in selected {
        let started = Instant::now();
        let topology = master_topology(loops, &master)?;
        let omega = topology.superficial_degree();
        let graph = topology.graph(&model, &mass)?;
        let mut record = serde_json::json!({
            "master": format!("{loops}L/{master}"),
            "loops": loops,
            "omega": omega,
            "edges": topology.edges,
            "extra_terms": extra,
        });
        let outcome = (|| -> Result<bool> {
            let t = Instant::now();
            let direct = massive_vacuum_direct(&graph, &model, &vakint, &settings, extra)?
                .series(GS.dim_epsilon, 0, -1)?
                .to_atom()
                .expand();
            record["direct_ms"] = (t.elapsed().as_millis() as u64).into();
            record["direct"] = direct.to_string().into();
            record["direct_unit"] = unit_scales(&direct, mass_symbol).to_string().into();
            let t = Instant::now();
            let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
            record["plan_ms"] = (t.elapsed().as_millis() as u64).into();
            record["spinneys"] = plan
                .ordinary_uv
                .iter()
                .filter(|term| !term.counterterms.is_empty())
                .count()
                .into();
            let sub = plan.massive_vacuum_subdivergence_poles(
                &graph, &model, &mass, &vakint, &settings, extra,
            )?;
            record["children_ms"] = (sub.children_ms as u64).into();
            record["subdivergences"] = sub.poles.to_string().into();
            record["terms"] = sub
                .terms
                .iter()
                .map(|term| {
                    serde_json::json!({
                        "label": term.label,
                        "child_loops": term.child_loops,
                        "quotient_loops": term.quotient_loops,
                        "ms": term.elapsed_ms as u64,
                        "contribution_poles": term.contribution
                            .series(GS.dim_epsilon, 0, -1)
                            .map(|s| s.to_atom().expand().to_string())
                            .unwrap_or_default(),
                    })
                })
                .collect::<Vec<_>>()
                .into();
            record["m_uv_in_subdivergences"] =
                sub.poles.contains_symbol(GS.m_uv_vacuum).into();
            let mut ok = true;
            let mut deltas = Vec::new();
            if routes.split(',').any(|r| r.trim() == "i") {
                let t = Instant::now();
                let delta = plan.massive_vacuum_counterterm(&graph, &model, &vakint, &settings)?;
                record["delta_i_ms"] = (t.elapsed().as_millis() as u64).into();
                record["delta_i"] = delta.to_string().into();
                record["m_uv_in_delta_i"] = delta.contains_symbol(GS.m_uv_vacuum).into();
                deltas.push(("i", delta));
            }
            if routes.split(',').any(|r| r.trim() == "ii") {
                let t = Instant::now();
                let (delta, log) =
                    topology.overall_counterterm_rstar(&model, &mass, &vakint, &settings)?;
                record["delta_ii_ms"] = (t.elapsed().as_millis() as u64).into();
                record["delta_ii"] = delta.to_string().into();
                record["delta_ii_terms"] = log
                    .iter()
                    .map(|(a, w, d, ms)| {
                        serde_json::json!({"a": a, "weight": w, "delta": d.to_string(), "ms": *ms as u64})
                    })
                    .collect::<Vec<_>>()
                    .into();
                deltas.push(("ii", delta));
            }
            for (route, delta) in &deltas {
                let identity = (-delta + &sub.poles).expand();
                let difference = (&identity - &direct).expand().together();
                let difference_unit = unit_scales(&difference, mass_symbol).together();
                record[format!("identity_{route}")] = identity.to_string().into();
                record[format!("difference_{route}")] = difference.to_string().into();
                record[format!("difference_unit_{route}")] = difference_unit.to_string().into();
                let pass = difference.is_zero();
                record[format!("pass_{route}")] = pass.into();
                ok &= pass;
            }
            if deltas.len() == 2 {
                record["delta_i_minus_ii"] =
                    (&deltas[0].1 - &deltas[1].1).expand().together().to_string().into();
            }
            Ok(ok)
        })();
        record["total_ms"] = (started.elapsed().as_millis() as u64).into();
        match outcome {
            Ok(ok) => {
                record["ok"] = ok.into();
                if !ok {
                    failures.push(format!("{loops}L/{master}"));
                }
            }
            Err(error) => {
                record["ok"] = false.into();
                record["error"] = format!("{error:?}").into();
                failures.push(format!("{loops}L/{master} (error)"));
            }
        }
        println!("RSTAR_MASSIVE {record}");
        if let Some(out) = &out {
            let mut file = std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(out)?;
            writeln!(file, "{record}")?;
        }
    }
    assert!(failures.is_empty(), "failed: {failures:?}");
    Ok(())
}
