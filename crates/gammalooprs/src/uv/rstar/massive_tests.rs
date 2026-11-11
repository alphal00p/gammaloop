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
        numerators: Vec::new(),
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
            record["children"] = sub
                .children
                .iter()
                .map(|(label, loops, edges, value)| {
                    serde_json::json!({"label": label, "loops": loops, "edges": edges, "m_uv": value.expand().contains_symbol(GS.m_uv_vacuum), "delta": value.expand().to_string()})
                })
                .collect::<Vec<_>>()
                .into();
            let implied = (&sub.poles - &direct).expand();
            record["delta_implied"] = implied.to_string().into();
            record["delta_implied_has_log"] =
                implied.contains_symbol(symbolica::atom::Symbol::LOG).into();
            record["delta_implied_has_m_uv"] = implied.contains_symbol(GS.m_uv_vacuum).into();
            record["m_uv_in_subdivergences"] = sub.poles.contains_symbol(GS.m_uv_vacuum).into();
            let mut ok = true;
            let mut deltas = Vec::new();
            if routes.split(',').any(|r| r.trim() == "i") {
                let t = Instant::now();
                match plan.massive_vacuum_counterterm(&graph, &model, &vakint, &settings) {
                    Ok(delta) => {
                        record["delta_i_ms"] = (t.elapsed().as_millis() as u64).into();
                        record["delta_i"] = delta.to_string().into();
                        record["m_uv_in_delta_i"] = delta.contains_symbol(GS.m_uv_vacuum).into();
                        deltas.push(("i", delta));
                    }
                    Err(error) => {
                        record["error_i"] = format!("{error:?}").into();
                        ok = false;
                    }
                }
            }
            if routes.split(',').any(|r| r.trim() == "ii") {
                let t = Instant::now();
                let (delta, log) =
                    match topology.overall_counterterm_rstar(&model, &mass, &vakint, &settings) {
                        Ok(value) => value,
                        Err(error) => {
                            record["error_ii"] = format!("{error:?}").into();
                            return Ok(false);
                        }
                    };
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
            for a in 0..deltas.len() {
                for b in a + 1..deltas.len() {
                    record[format!("delta_{}_minus_{}", deltas[a].0, deltas[b].0)] = (&deltas[a].1
                        - &deltas[b].1)
                        .expand()
                        .together()
                        .to_string()
                        .into();
                }
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

/// Part (A) only, `-K sum_S Delta(S) G/S`, for topologies listed in the JSON
/// file `RSTAR_PARTA_SPEC` (`[{"name":..,"loops":..,"edges":[[a,b,p],..]}]`,
/// scalar numerators only), optionally restricted by `RSTAR_PARTA_NAMES`.
/// Records (JSON lines) go to `RSTAR_PARTA_OUT`; each holds the Laurent
/// coefficients eps^-L..eps^-1 of (A) at m = mu = 1 and per-term timings.
#[test]
#[ignore = "needs VAKINT_TEST_EXTERNAL_DATA and RSTAR_PARTA_SPEC"]
fn massive_vacuum_part_a_from_spec() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    let settings = super::tests::settings();
    let spec: serde_json::Value = serde_json::from_slice(&std::fs::read(
        std::env::var("RSTAR_PARTA_SPEC").expect("RSTAR_PARTA_SPEC"),
    )?)?;
    let names = std::env::var("RSTAR_PARTA_NAMES").ok().map(|v| {
        v.split(',')
            .map(|s| s.trim().to_owned())
            .collect::<Vec<_>>()
    });
    let extra: i64 = std::env::var("RSTAR_MASSIVE_EXTRA")
        .ok()
        .map(|v| v.parse().unwrap())
        .unwrap_or(0);
    let out = PathBuf::from(std::env::var("RSTAR_PARTA_OUT").expect("RSTAR_PARTA_OUT"));
    let vakint = loaded_vakint(4)?;
    let mass_symbol = symbol!("rstar_vacuum_mass"; Scalar);
    let mass = Atom::var(mass_symbol);
    let mut failures = Vec::new();
    for entry in spec.as_array().expect("a JSON list") {
        let name = entry["name"].as_str().unwrap().to_owned();
        if names.as_ref().is_some_and(|names| !names.contains(&name)) {
            continue;
        }
        let topology = VacuumTopology {
            name: format!("m5L_{name}"),
            edges: entry["edges"]
                .as_array()
                .unwrap()
                .iter()
                .map(|e| {
                    (
                        e[0].as_u64().unwrap() as usize,
                        e[1].as_u64().unwrap() as usize,
                        e[2].as_u64().unwrap() as u32,
                    )
                })
                .collect(),
            numerators: entry["numerators"]
                .as_array()
                .map(|factors| {
                    factors
                        .iter()
                        .map(|factor| {
                            factor
                                .as_array()
                                .unwrap()
                                .iter()
                                .map(|t| (t[0].as_u64().unwrap() as usize, t[1].as_i64().unwrap()))
                                .collect()
                        })
                        .collect()
                })
                .unwrap_or_default(),
        };
        let started = Instant::now();
        let mut record = serde_json::json!({
            "numerators": topology.numerators,
            "master": name,
            "loops": topology.loops(),
            "omega": topology.superficial_degree(),
            "edges": topology.edges,
            "extra_terms": extra,
        });
        let outcome = (|| -> Result<()> {
            let monomials = topology.numerator_monomials(&mass);
            record["monomials"] = monomials.len().into();
            let validate = std::env::var("RSTAR_PARTA_VALIDATE").is_ok_and(|v| v == "1");
            let mut poles = Atom::Zero;
            let mut direct = Atom::Zero;
            let mut delta_i = Atom::Zero;
            let mut spinneys = 0usize;
            let mut children = Vec::new();
            let mut terms = Vec::new();
            let mut children_ms = 0u128;
            for monomial in &monomials {
                let graph = topology.monomial_graph(&model, &mass, &monomial.attachments)?;
                let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
                spinneys += plan
                    .ordinary_uv
                    .iter()
                    .filter(|term| !term.counterterms.is_empty())
                    .count();
                let sub = plan.massive_vacuum_subdivergence_poles(
                    &graph, &model, &mass, &vakint, &settings, extra,
                )?;
                poles += &monomial.coefficient * &sub.poles;
                children.extend(sub.children);
                terms.extend(sub.terms);
                children_ms += sub.children_ms;
                if validate {
                    direct += &monomial.coefficient
                        * massive_vacuum_direct(&graph, &model, &vakint, &settings, extra)?;
                    delta_i += &monomial.coefficient
                        * plan.massive_vacuum_counterterm(&graph, &model, &vakint, &settings)?;
                }
            }
            record["spinneys"] = spinneys.into();
            let sub = super::massive::MassiveVacuumSubdivergences {
                poles: poles.series(GS.dim_epsilon, 0, -1)?.to_atom().expand(),
                terms,
                children_ms,
                children,
            };
            if validate {
                let direct = direct.series(GS.dim_epsilon, 0, -1)?.to_atom().expand();
                let delta_i = delta_i.expand();
                let difference = (&sub.poles - &delta_i - &direct).expand().together();
                record["direct"] = direct.to_string().into();
                record["delta_i"] = delta_i.to_string().into();
                record["m_uv_in_delta_i"] = delta_i.contains_symbol(GS.m_uv_vacuum).into();
                record["difference_i"] = difference.to_string().into();
                record["pass_i"] = difference.is_zero().into();
            }
            record["children_ms"] = (sub.children_ms as u64).into();
            record["part_a"] = sub.poles.to_string().into();
            let unit = unit_scales(&sub.poles, mass_symbol);
            record["part_a_unit"] = unit.to_string().into();
            let series = unit.series(GS.dim_epsilon, 0, -1)?;
            let mut coefficients = serde_json::Map::new();
            let mut exact_file = serde_json::Map::new();
            for (power, coefficient) in series.terms() {
                coefficients.insert(power.to_string(), coefficient.expand().to_string().into());
                exact_file.insert(
                    power.to_string(),
                    symbolica::printer::AtomPrinter::new_with_options(
                        coefficient.expand().as_view(),
                        symbolica::printer::PrintOptions::file(),
                    )
                    .to_string()
                    .into(),
                );
            }
            record["part_a_coefficients"] = coefficients.into();
            record["part_a_coefficients_file"] = exact_file.into();
            let mut numeric_settings = settings.vakint.true_settings();
            let digits: u32 = std::env::var("RSTAR_NUMERIC_DIGITS")
                .ok()
                .map(|v| v.parse().unwrap())
                .unwrap_or(120);
            numeric_settings.run_time_decimal_precision = digits;
            record["part_a_numeric_digits"] = digits.into();
            match vakint.numerical_evaluation(
                &numeric_settings,
                unit.as_view(),
                &Default::default(),
                &Default::default(),
                None,
            ) {
                Ok((value, _)) => {
                    let mut numeric = serde_json::Map::new();
                    for (power, coefficient) in value.coefficients() {
                        numeric.insert(
                            power.to_string(),
                            format!("{:.*e}", digits as usize + 5, coefficient.value().re).into(),
                        );
                        numeric.insert(
                            format!("{power}_im"),
                            format!("{:.5e}", coefficient.value().im).into(),
                        );
                    }
                    record["part_a_numeric"] = numeric.into();
                }
                Err(error) => record["part_a_numeric_error"] = format!("{error}").into(),
            }
            record["m_uv_in_part_a"] = sub.poles.contains_symbol(GS.m_uv_vacuum).into();
            if std::env::var("RSTAR_PARTA_DELTA").is_ok_and(|v| v == "ii") {
                let t = Instant::now();
                match topology.overall_counterterm_rstar(&model, &mass, &vakint, &settings) {
                    Ok((delta, log)) => {
                        let unit = unit_scales(&delta, mass_symbol);
                        record["delta_ii"] = unit.to_string().into();
                        let mut coefficients = serde_json::Map::new();
                        let mut numeric = serde_json::Map::new();
                        for (power, coefficient) in unit.series(GS.dim_epsilon, 0, -1)?.terms() {
                            coefficients
                                .insert(power.to_string(), coefficient.expand().to_string().into());
                        }
                        let mut numeric_settings = settings.vakint.true_settings();
                        numeric_settings.run_time_decimal_precision = 120;
                        if let Ok((value, _)) = vakint.numerical_evaluation(
                            &numeric_settings,
                            unit.as_view(),
                            &Default::default(),
                            &Default::default(),
                            None,
                        ) {
                            for (power, coefficient) in value.coefficients() {
                                numeric.insert(
                                    power.to_string(),
                                    format!("{:.125e}", coefficient.value().re).into(),
                                );
                            }
                        }
                        record["delta_ii_coefficients"] = coefficients.into();
                        record["delta_ii_numeric"] = numeric.into();
                        record["delta_ii_graphs"] = log.len().into();
                        if validate {
                            let direct = direct.series(GS.dim_epsilon, 0, -1)?.to_atom().expand();
                            let difference = (&sub.poles - &delta - &direct).expand().together();
                            record["difference_ii"] = difference.to_string().into();
                            record["pass_ii"] = difference.is_zero().into();
                        }
                        record["delta_ii_ms"] = (t.elapsed().as_millis() as u64).into();
                        record["delta_ii_terms"] = log
                            .iter()
                            .map(|(a, w, d, ms)| serde_json::json!({"a": a, "weight": w, "delta": d.to_string(), "ms": *ms as u64}))
                            .collect::<Vec<_>>()
                            .into();
                    }
                    Err(error) => record["delta_ii_error"] = format!("{error:?}").into(),
                }
            }
            // Child polynomials are written in all boundary momentum labels without
            // imposing momentum conservation, so terms that vanish by conservation
            // may carry m_uv; the quotient integration imposes it. Check every
            // integrated term contribution (and the total) instead.
            record["terms_with_m_uv"] = sub
                .terms
                .iter()
                .filter(|term| term.contribution.expand().contains_symbol(GS.m_uv_vacuum))
                .map(|term| term.label.clone())
                .collect::<Vec<_>>()
                .into();
            record["terms"] = sub
                .terms
                .iter()
                .map(|term| {
                    serde_json::json!({
                        "label": term.label,
                        "child_loops": term.child_loops,
                        "quotient_loops": term.quotient_loops,
                        "ms": term.elapsed_ms as u64,
                    })
                })
                .collect::<Vec<_>>()
                .into();
            Ok(())
        })();
        record["total_ms"] = (started.elapsed().as_millis() as u64).into();
        let mismatch = record["pass_i"].as_bool() == Some(false)
            || record["pass_ii"].as_bool() == Some(false)
            || record["m_uv_in_part_a"].as_bool() == Some(true)
            || record["terms_with_m_uv"]
                .as_array()
                .is_some_and(|terms| !terms.is_empty())
            || record.get("delta_ii_error").is_some()
            || record.get("part_a_numeric_error").is_some();
        record["ok"] = (outcome.is_ok() && !mismatch).into();
        if let Err(error) = outcome {
            record["error"] = format!("{error:?}").into();
        }
        if record["ok"] != true {
            failures.push(name.clone());
        }
        println!("RSTAR_PARTA {record}");
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&out)?;
        writeln!(file, "{record}")?;
    }
    assert!(
        failures.is_empty(),
        "part (A) mismatches or errors: {failures:?}"
    );
    Ok(())
}
